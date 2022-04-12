module response
  use EOM_IMSRG
  use basic_IMSRG
  use cross_coupled
  implicit none

contains
!==============================================================================================
!==============================================================================================
subroutine compute_response_function(jbas,HS,OP,OPhu)
  
  integer :: N 
  integer :: nev
  type(spd) :: jbas
  type(sq_op) :: op,V1,Q1,Q2,w1,w2,PIVOT,HS,OPhu
  type(cc_mat) :: OpCC,QCC,WCC
  type(ex_pandya_mat) :: QPP,WPP
  type(ex_cc_mat) :: OpPP
  real(8),allocatable,dimension(:) :: workl,D,resid,workD
  real(8),allocatable,dimension(:,:) :: V,Z,VX,XXX 
  integer :: i,j,k,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2,strength,sm
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct
  real(8) :: temp1,temp2,temp3,fac1,fac2,RSD,RSD_prev
  integer :: iprev,iparam_c(11),ipntr_c(11),info_c,ido_c
  real(8),allocatable,dimension(:) :: workl_c,resid_c,workD_c
  real(8),allocatable,dimension(:,:) :: V_c

  fac1 = 1.39626338*197.326968*2.0/137.03599976
  fac2 = 1.39626338*3.1415926*3.1415926*4.0/137.03599976
  RSD = 0.d0
  RSD_prev = 0.d0

  open (unit=82, file = trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_gaussianfold_'//&
       OP%trans_label//'_response.dat')
  
  open (unit=83, file = trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_discrete_'//&
       OP%trans_label//'_response.dat')
    
!!! this is how many valid eigenvalues you seek, multiply by ten for
!!! number of pivots 
!  nev = 150
  nev = 500
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Q1%pphh_ph=.true.
  Q2%pphh_ph=.true.
  w1%pphh_ph=.false.
  w2%pphh_ph=.false.
  
  call duplicate_sq_op(OPhu,w1,'w') !workspace
  call duplicate_sq_op(OPhu,w2,'w') !workspace
  call duplicate_sq_op(OPhu,Q1,'y') !workspace
  call duplicate_sq_op(OPhu,Q2,'y') !workspace
  
  call init_ph_mat(Q1,OpPP,jbas) !cross coupled ME
  call init_ph_mat(Q1,QPP,jbas) !cross coupled ME
  call init_ph_wkspc(QPP,WPP) 
  
  h = HS%belowEF !holes
  p = HS%Nsp-h  !particles

  print*
  print*
  print*, "==========================================="
  print*, " COMPUTING "//OP%trans_label//" RESPONSE FUNCTION FOR "&
       //nucleus_name(HS%Aneut,HS%Aprot)
  print*, "==========================================="
  
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
  
        if (triangle(jbas%jj(i),jbas%jj(j),OP%rank)) then  

           if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
           if (mod(jbas%ll(i) + jbas%ll(j) + OP%dpar/2,2) .ne. 0 ) cycle
                     
           sps = sps+1
        end if 
     end do
  end do
  

  ! tensor case
  tps = 0 
  do q = 1, OP%nblocks
     
     do II = 1,size( OP%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( OP%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(OP%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   OP%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if
           
           if (mod(OP%tblck(q)%Jpair(2)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   OP%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if
           
           tps = tps+ 1
        end do
     end do
     
     if (OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2)) cycle
     
     do II = 1,size( OP%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( OP%tblck(q)%tgam(7)%X(1,:) )  
           
           if (mod(OP%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   OP%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if
           
           if (mod(OP%tblck(q)%Jpair(2)/2,2) == 1) then 
              
              if ( OP%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   OP%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if
           
           tps = tps+ 1
        end do
     end do
     
  end do

!bhu TDA
  if((sps+tps)>nTDA) tps=0
!bhu TDA
           
  print*, '1p1h Amplitudes: ', sps
  print*, '2p2h Amplitudes: ', tps
  N = sps + tps ! number of ph and pphh SDs 
  Q1%neq = N
  Q2%neq = N
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SM' ! compute smallest eigenvalues in magnitude ('SA' is algebraic). 
  tol = 1.0E-6 ! error tolerance?
  info = 1 ! THIS MEANS USE THE PIVOT I TELL YOU TO. 
!  info = 0
  ncv = nev+1 ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),VX(N,2),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 1   !!!! THIS TURNS OFF IMPLICIT RESTART
!bhu  mxiter = 300
  mode = 1
   
  allocate(resid(N),workD(3*N)) 
  allocate(resid_c(N),workD_c(3*N)) 
  allocate(V_c(N,NCV),workl_c(lworkl))

  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  call rewrap_tensor_pphh( V(:,1), OP  ,N ,jbas) !!! GET PIVOT

  ! normalize pivot
  strength = sum(V(:,1)**2)
  V(:,1) = V(:,1)/sqrt(strength) 
  resid = V(:,1)

  VX(:,1) = V(:,1) 

  rvec= .true. 
  howmny = 'A'
  allocate(selct(NCV)) 
  allocate(D(NEV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  iparam_c(5) = -1

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this. 
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
          ncv, v, ldv, iparam, ipntr, workd, workl, &
          lworkl, info )
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )

     write(*,*) ', NCONV: ',iparam(5)
     if(i<NEV) then
        write(82,'(i6,a,i4)',advance='yes') i, ', NCONV: ',iparam(5)
        BACKSPACE(UNIT=82)
     end if

     resid_c = resid
     v_c = v
     iparam_c = iparam
     ipntr_c = ipntr
     workd_c = workd
     workl_c = workl
     info_c = info

     iparam(5) = i
     if(i>NEV) iparam(5) = NEV
     
     if(i >= nev-1) then    
     call dseupd( rvec, howmny, selct, d, Z, ldz, sigma, &
          bmat, n, which, nev, tol, resid, ncv, v, ldv, &
          iparam, ipntr, workd, workl, lworkl, info )  


     if(i==NEV) then
        do j = 1, NEV
           temp1 = strength * sum(Z(:,j)*VX(:,1))**2
           write(*,'(I3,a,2(f16.8))') j,' Ev, RSG:', D(j), temp1
           write(83,'(2(f25.14))') D(j), temp1
        end do
     end if
     RSD = 0.d0
     do j = 1, NEV
        RSD = RSD + fac1 * strength * sum(Z(:,j)*VX(:,1))**2 / D(j)
     end do
     write(*,*) 'response strength is: ', RSD
     write(83,*) i, iparam_c(5), ' response strength is: ', RSD
     if(abs(RSD-RSD_prev)<1.0E-6 .and. i > 2*NCV) cycle
     RSD_prev = RSD
     end if
     
     resid = resid_c
     v = v_c
     iparam = iparam_c
     ipntr = ipntr_c
     workd = workd_c
     workl = workl_c
     info = info_c

     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if

     if ( i > 10*NCV ) STOP 'response failed to converge' 
     call matvec_nonzeroX_prod(N,HS,Q1,Q2,w1,w2,OpPP,QPP,WPP,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 

     i=i+1 

  end do
  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )  

!  x = 0.d0 
!  do i = 1, 10000
!     sm = 0.d0 
!     do j = 1,NEV
!        sm = sm + strength * sum(Z(:,j)*VX(:,1))**2 * gaussian(x,D(j),0.5d0)        
!     end do     
!     write(82,'(2(f25.14))') x ,  sm
!     x = x + .01d0
!  end do

  temp1=0.d0
  temp2=0.d0
  temp3=0.d0
  RSD_prev=0.d0
  do j = 1, NEV
    RSD = strength * sum(Z(:,j)*VX(:,1))**2
    temp1 = temp1 + D(j) * RSD
    temp2 = temp2 + RSD
    temp3 = temp3 + fac1 * RSD / D(j)
    write(*,'(i3,2(f16.8))') j, D(j), RSD
    write(83,'(2(f25.14))') D(j), RSD
!bhu    call rewrap_tensor_pphh( VX(:,2), OP  ,N ,jbas)
!bhu    write(*,'(i3,3(f16.8))') j, D(j),  sum(Z(:,j)*VX(:,2))**2,&
!bhu    strength * sum(Z(:,j)*VX(:,1))**2
  end do

  write(*,*) 'electric dipole polarizability: ', temp3
  write(83,*) 'electric dipole polarizability: ', temp3
  write(83,*) 'm1, energy weight E*R sum:      ', temp1
  write(83,*) 'm0, R(E) sum:                   ', temp2
  write(83,*) 'centroid m1/m0:                 ', temp1/temp2
  write(83,'(2(a20),a40)')'E    ','R(E)    ','photodisintegration cross section'

  close(82)
  close(83)
  
end subroutine COMPUTE_RESPONSE_FUNCTION

real(8) function gaussian(x,E,sig)
  implicit none

  real(8),intent(in) :: E,sig,x
  real(8) :: div
  
  div = sqrt(2*PI_const)*sig 
  gaussian = exp( -(x-E)**2 /(2.d0*sig**2) )/div 
end function gaussian

real(8) function dirac_delta(x,E,sig)
  implicit none

  real(8),intent(in) :: E,sig,x
  
  if ( abs(x-E) >1e-3) then 
     dirac_delta = 0.0
  else
     dirac_delta = 1.0
  end if
  
end function dirac_delta


subroutine rewrap_tensor_pphh( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  real(8) :: pre

  pre = 1.d0
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  
  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           v(i) = AX%fph(IX,JX)
           i = i + 1
        end if 
     end do
  end do
  
!bhu TDA
  if(i>nTDA) return
!bhu TDA

  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,size( AX%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%tgam(3)%X(II,JJ)

           if( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
               AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) then
              v(i) = v(i) * pre
           end if
           if( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
               AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) then
              v(i) = v(i) * pre
           end if

           i = i + 1
        end do
     end do

     ! IF THE Js ARE THE SAME THEN WE NEED TO MAKE SURE THAT IS REFLECTED IN 
     ! THE TRANSPOSE
     if (AX%tblck(q)%Jpair(1) == AX%tblck(q)%Jpair(2)) then 
        
        cycle 
     
     end if
     
     ! OTHERWISE BUSINESS AS USUAL.
     do II = 1,size( AX%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if

           v(i)=AX%tblck(q)%tgam(7)%X(II,JJ)*AX%tblck(q)%lam(1)*AX%herm

           if( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
               AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) then
              v(i) = v(i) * pre
           end if
           if( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
               AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) then
              v(i) = v(i) * pre
           end if

           i = i + 1
        end do
     end do

  end do

end subroutine rewrap_tensor_pphh



end module response


