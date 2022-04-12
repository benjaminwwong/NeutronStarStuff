program main_IMSRG
!!!===================================================================
!!!     Equations of Motion for closed-shell nuclei 
!!!====================================================================
  use isospin_operators
  use HF_mod
  use response
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  use interaction_IO
  use EOM_IMSRG
  use brute_force_testing
  use three_body_routines
  use deuteron
  implicit none
  
  type(spd) :: jbas,jbx
  type(sq_op) :: HS,ETA,DH,w1,w2,rirj,pipj,r2_rms,Otrans,exp_omega
  type(sq_op) :: num,cr,H0,Hcm,Oscal 
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  type(iso_ladder),allocatable,dimension(:) :: isoladder_ops
  type(iso_operator) :: GT_trans
  type(cc_mat) :: CCHS,CCETA,WCC
  type(full_sp_block_mat) :: coefs,coefsT,TDA,ppTDA,rrTDA
  type(three_body_force) :: threebod
  type(obsv_mgr) :: trans,moments
  type(eom_mgr) :: eom_states
  character(200) :: inputs_from_command
  character(10) :: other_obs
  character(1) :: quads,trips,trans_type
  integer :: i,j,T,JTot,a,b,c,d,g,q,ham_type,j3,ix,jx,kx,lx,PAR,Tz,trans_rank
  integer :: np,nh,nb,k,l,m,n,method_int,mi,mj,ma,mb,j_min,ex_Calc_int,J1,J2
  integer :: na,la,lb,totstates,numstates,oldnum,qx,dTZ,oldnum_dTz,numstates_dTz
  real(8) :: sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx,dcgi,dcgi00
  real(8) :: corr_op,pre,x,corr,de_trips
  logical :: hartree_fock,COM_calc,r2rms_calc,me2j,me2b,trans_calc
  logical :: skip_setup,skip_gs,do_HF,TEST_commutators,mortbin,decouple
  external :: build_gs_white,build_specific_space,build_gs_atan,build_gs_w2
  external :: build_ex_imtime,build_sd_shellmodel

!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  time0 = omp_get_wtime()

  call getarg(1,inputs_from_command) 
  call getarg(2,resubmitter) 
  
  if (trim(inputs_from_command) == 'X') then 
     test_commutators = .true.
     inputs_from_command = ''
  else
     test_commutators = .false.
  end if

  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,ex_calc_int,COM_calc,r2rms_calc,other_obs,me2j,&
       me2b,mortbin,jbas%hw,skip_setup,skip_gs,quads,trips,threebod%e3max)

  call read_tokyo_orbits(jbas,HS%Aprot,HS%Aneut,trips,jbx)

  call print_system(jbas) 

  call allocate_blocks(jbas,HS)

  HS%herm = 1
  HS%hospace = jbas%hw
  
  do_hf = .true. 

!=============================================================
! READ HAMILTONIAN
!=============================================================
  print*, 'reading 2-body Hamiltonian'
  call read_tokyo_hamiltonian(HS,jbas)

  call print_time

  call print_header
     
  decouple=.false.


!=======================================================================
!  equations-of-motion calculation 
!=======================================================================
91 if (ex_calc_int==1) then
     print*, 'STARTING EXCITED STATES CALCULATION'
     totstates=read_eom_file(trans,moments,eom_states,jbas)! total number of states
     allocate(ladder_ops(totstates-eom_states%total_dTz))
     allocate(isoladder_ops(eom_states%total_dTz))     
     
     oldnum = 0
     oldnum_dTz = 0
     numstates = 0
     numstates_dTz = 0
     do q = 1,eom_states%num

        if (eom_states%dTz(q) == 0 ) then 
           oldnum = oldnum + Numstates
           Numstates = eom_states%number_requested(q)        
           ladder_ops(1+oldnum:Numstates+oldnum)%xindx = q
           call calculate_excited_states(eom_states%ang_mom(q),eom_states%par(q),numstates,HS,&
                jbas,ladder_ops(1+oldnum:Numstates+oldnum))
        
           call print_time
           if (eom_states%trips) then
              call print_triples_header
              do qx = 1+oldnum,Numstates+oldnum
                 
                 t1= omp_get_wtime()
                 dE_trips=EOM_triples(HS,ladder_ops(qx),jbas)  
                 t2= omp_get_wtime()
                                  
                 write(*,'(A2,4(f20.10))') eom_states%name(q),&
                      ladder_ops(qx)%E0,ladder_ops(qx)%E0+dE_trips,dE_trips,t2-t1
              end do
           end if

        else
   
           oldnum_dTz = oldnum_dTz + Numstates_dTz
           Numstates_dTz = eom_states%number_requested(q)        
           isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz)%xindx = q
           call calculate_isospin_states(eom_states%ang_mom(q),eom_states%par(q),eom_states%dTz(q),&
                numstates_dTZ,HS,jbas,isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz))
        
        end if        
               
    end do

    if ((trans%num + moments%num)*(Numstates+Numstates_dTz) == 0) then 
       trans_Calc=.false.
    else
       trans_Calc= .true.
    end if
    
    call print_time
    
     Otrans%xindx = eom_states%num+1
     GT_Trans%xindx = Otrans%xindx

     trans_type = trans%oper(1:1)

     if (trans_type == 'G') then 
        call initialize_transition_operator('G',trans%oper(2:2),GT_trans,HS,jbas,trans_calc)
     else        
        read(trans%oper(2:2),'(I1)') trans_rank        
        call initialize_transition_operator(trans_type,trans_rank,Otrans,HS,jbas,trans_calc)
     end if 

     if (writing_human) then
        call write_onebody_tensor_human(Otrans,jbas,"E2_0")   
        call write_twobody_tensor_human(Otrans,jbas,"E2_0")
     end if

     if (trans_calc) then 
        print*, 'reading transition operator...' 

        if (trans_type == 'G') then 
           print*, 'iso_operator is not implemented!'                                 
           call read_operator_human(Otrans,jbas)                                 
        else
           call read_operator_human(Otrans,jbas)                                 
        end if
   
     end if
     call print_time

     
     if (trans_type == 'G') then 
        call EOM_beta_observables( ladder_ops, isoladder_ops, GT_trans, HS, Hcm,trans, moments,eom_states,jbas)
     else           
        call EOM_observables( ladder_ops, isoladder_ops, Otrans, HS, Hcm,trans, moments,eom_states,jbas)

        if (eom_states%response) then 
!bhu           call compute_response_function(jbas,HS,Otrans) 
!bhu           call compute_response_function(jbas,HS,ladder_ops(1)) 
           call compute_response_function(jbas,HS,Otrans,ladder_ops(1)) 
        end if
        
     end if
     deallocate(isoladder_ops,ladder_ops)
  end if

  call print_time
contains

subroutine test
  ! call compare_tensor_scalar_commutator(jbas,-1,1) 
  ! stop
  ! deallocate(jbas%xmap) 
  ! call test_scalar_scalar_commutator(jbas,-1,1) 
  ! deallocate(jbas%xmap)
  ! call test_EOM_scalar_scalar_commutator(jbas,1,1)
  ! deallocate(jbas%xmap)
!   call test_EOM_scalar_tensor_commutator(jbas,1,1,6,2)  
!   deallocate(jbas%xmap,jbas%xmap_tensor,phase_hh,phase_pp)
!   deallocate(half6j%tp_mat)
!  call test_scalar_tensor_commutator(jbas,-1,1,6,2) 
 ! call test_tensor_product(jbas,1,1,2,2,2,2,0,2) 
!  call test_EOM_iso_commutator(jbas,1,1,4,0,0)
!  call test_scalar_iso_commutator(jbas,-1,1,6,2,1) 
  call test_tensor_dTZ_product(jbas,1,1,2,2,2,2,0,2,-1) 

end subroutine test
end program main_IMSRG
!=========================================================================
subroutine print_header
  implicit none 
  
  print*, '================================'//&
       '==================================='
  print*, '  iter        s            E0      '//&
       '    E0+MBPT(2)      |MBPT(2)|  '  
  print*, '================================'//&
       '==================================='

end subroutine print_header

subroutine print_triples_header
  implicit none
  
  print*
  print*, '====================================='//&
       '=============================================='
  print*, 'J^Pi      dE(1+2)            dE(1+2+3)'//&
       '             dE(3)               time    '
  print*, '====================================='//&
       '=============================================='
  
end subroutine print_triples_header
