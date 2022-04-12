#!/usr/bin/env python

##########################################################################
##  goUniversal.py
##
##  A python script to run or submit jobs for the common use cases
##  of the IMSRG++ code. We check whether there is a pbs or slurm
##  scheduler, assign the relevant input parameters, set names
##  for the output files, and run or submit.
##  						-Ragnar Stroberg
##  						TRIUMF Nov 2016
######################################################################

from os import path,environ,mkdir,remove
from sys import argv
from subprocess import call,PIPE
from time import time,sleep
from datetime import datetime

### Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'SLURM'

### The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS=32
exe = '%s/bin/imsrg++'%(environ['HOME'])

### Flag to swith between submitting to the scheduler or running in the current shell
#batch_mode=False
batch_mode=True
if 'terminal' in argv[1:]: batch_mode=False

### Don't forget to change this. I don't want emails about your calculations...
mail_address = 'bhu@triumf.ca'

### This comes in handy if you want to loop over Z
ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}

### Maximum value of s, and maximum step size ds
ARGS['smax'] = '500'
ARGS['dsmax'] = '0.5'

#ARGS['lmax3'] = '10' # for comparing with Heiko

### Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.25'

### Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
#ARGS['scratch'] = 'SCRATCH'    
#ARGS['scratch'] = '/global/scratch/bhu/scratch'    

### Generator for core decoupling, can be atan, white, imaginary-time.  (atan is default)
#ARGS['core_generator'] = 'imaginary-time' 
### Generator for valence deoupling, can be shell-model, shell-model-atan, shell-model-npnh, shell-model-imaginary-time (shell-model-atan is default)
#ARGS['valence_generator'] = 'shell-model-imaginary-time' 

### Solution method
ARGS['method'] = 'magnus'
#ARGS['method'] = 'NSmagnus'
#ARGS['method'] = 'brueckner'
#ARGS['method'] = 'flow'
#ARGS['method'] = 'HF'
#ARGS['method'] = 'MP3'

#ARGS['write_omega'] = 'true'

### Tolerance for ODE solver if using flow solution method
#ARGS['ode_tolerance'] = '1e-5'

if BATCHSYS == 'PBS':
  FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q oak
#PBS -d %s
#PBS -l walltime=288:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=250gb
#PBS -m ae
#PBS -M %s
#PBS -j oe
#PBS -o imsrg_log/%s.o
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
%s
  """

elif BATCHSYS == 'SLURM':
  FILECONTENT = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=imsrg_log/%s.%%j
#SBATCH --time=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""

### Make a directory for the log files, if it doesn't already exist
if not path.exists('imsrg_log'): mkdir('imsrg_log')

#As=[19,23,27,29,73]
#Zs=[9,11,13,14,32]
As=[19,23,27,29]
Zs=[9,11,13,14]
#As=[16]
#Zs=[8]
#As=[19]
#Zs=[9]
#As=[23]
#Zs=[11]
#As=[27]
#Zs=[13]
#As=[29]
#Zs=[14]
#As=[73]
#Zs=[32]
#As=[127]
#Zs=[53]
#As=[129]
#Zs=[54]
#As=[131]
#Zs=[54]

### Loop over multiple jobs to submit
#for A in range(9,10):
#for Z in [9]:
 #A=19
for i in range(0,len(As)):
 Z=Zs[i]
 A=As[i]
 for reference in ['%s%d'%(ELEM[Z],A)]:
  ARGS['reference'] = reference
  print('Z = ',Z)
  #A=19
  for e in [4]:
   for hw in [20]:
#   for hw in [13.3530]:    #F19
#   for hw in [12.2222]:
#   for hw in [20]:

     ARGS['emax'] = '%d'%e
     ARGS['e3max'] = '0'
#     ARGS['emax_imsrg'] = '4'
#     ARGS['e2max_imsrg'] = '8'
#     ARGS['e3max_imsrg'] = '0'

#     ARGS['core_generator'] = 'imaginary-time'
#     ARGS['valence_generator'] = 'shell-model-imaginary-time'
#     ARGS['core_generator'] = 'white'
#     ARGS['valence_generator'] = 'shell-model'
#     ARGS['method'] = method

     ARGS['valence_file_format'] = 'tokyo'
#     ARGS['denominator_delta_orbit'] = 'all'
#     ARGS['denominator_delta'] = '10'
#     ARGS['goose_tank'] = 'true'

#     ARGS['IMSRG3'] = 'true'
#     ARGS['imsrg3_n7'] = 'true'
#     ARGS['dE3max'] = '10'
#     ARGS['OccNat3Cut'] = '0.00009'
#     ARGS['only_2b_eta'] = 'false'
#     ARGS['only_2b_omega'] = 'false'

#     ARGS['occ_file'] = 'occ_file_He9.dat'
#     ARGS['freeze_occupations'] = 'true'
#     ARGS['BetaCM'] = '5'
     ARGS['basis'] = 'NAT'
#     ARGS['use_NAT_occupations'] = 'true'
#     ARGS['basis'] = 'oscillator'
#     ARGS['hwBetaCM'] = '12'
#     ARGS['eta_criterion'] = '1e-4'

### Gamow IMSRG
#     inputpath = '/global/scratch/bhu/GHF/nucl-k/'
#     man='k-n-d3_pipj_n30_test2_nnlo-opt_N04_hw%d.0_%sA%d'%(hw,reference,A)
#     inputpath = '/global/scratch/bhu/GHF/nucl-MB/'
#     inputpath = '/global/scratch/bhu/GHF/output/'
#     man1 = 'n-p_n15_hagen_N3LO_srg2.15'
#     man1 = 'n-p_n24_hagen_N4LO_srg2.4'
#     man1 = 'p-sd_n24_zh_EM1.8_2.0'
#     man1 = 'n-d3f7p_n15_zh_N3LO_3N400'
#     man1 = 'N2LO_opt'
#     man1 = 'n-d3f7p_n15_test2_nnlo-opt'
#     man = man1 + '_%s_e%d_E%s_hw%d_A%d'%(reference,e,ARGS['e3max'],hw,A)
#     ARGS['inputorbits'] = inputpath+man+'_orb.dat'
#     ARGS['input1bme'] = inputpath+man+'_1B.dat'
#     ARGS['2bme'] = inputpath+man+'_2B.dat'
#     ARGS['3bme'] = 'none'
#     ARGS['LECs'] = man1+"_O22"
#     ARGS['basis'] = 'oscillator'
### GAmow IMSRG

#     ARGS['fmt2'] = 'oslo'
#     ARGS['LECs'] = 'N3LO'
#     ARGS['2bme'] = '/global/scratch/bhu/renorm/renorm_old/n3lo_N12_hw%d.dat'%(hw)
#     ARGS['3bme'] = 'none'

#     ARGS['LECs'] = 'gst_N3LO'
#     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
#     ARGS['2bme'] = '/global/scratch/exch/n3lo_int/TBMEA2n3lo-srg2.0_14_28.%d_TUD.int.gz'%(hw)
#     ARGS['fmt3'] = 'navratil'
#     ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
#     ARGS['3bme'] = '/global/scratch/exch/n3lo_int/v3trans_J3T3.int_NNn3lo_3Nind-srg2.0_220_161616.%d.gz'%(hw)
#     ARGS['3bme'] = '/global/scratch/exch/n3lo_int/v3trans_J3T3.int_NNn3lo_3NcD0.83cE-0.052_0.83-srg2.0_220_161616.%d.gz'%(hw)

#     ARGS['fmt2'] = 'oslo'
#     ARGS['LECs'] = 'N2LO_opt'
#     ARGS['2bme'] = '/global/scratch/bhu/renorm/renorm_nnlo/nnlo-opt_N12_hw%d.dat'%(hw)
#     ARGS['3bme'] = 'none'

#     ARGS['LECs'] = 'OccFile_NAT_N2LO_opt'
#     ARGS['2bme'] = '/global/scratch/exch/n2loopt_int/TBMEA2n2loopt_14_28.%d_TUD.int.gz'%(hw)
#     ARGS['3bme'] = 'none'

#     ARGS['LECs'] = 'EM1.8_2.0'
##     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
##     ARGS['2bme'] = '/global/scratch/exch/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
##     ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
##     ARGS['3bme'] = '/global/scratch/exch/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
#     ARGS['file2e1max'] = '18 file2e2max=36 file2lmax=18'
#     ARGS['2bme'] = '/global/scratch/exch/me2j/takayuki/TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw%d_emax18_e2max36.me2j.gz'%(hw)
#     ARGS['file3e1max'] = '18 file3e2max=36 file3e3max=24'
#     ARGS['3bme'] = '/global/scratch/exch/me3j/takayuki/NO2B_compact/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin'%(hw)
#     ARGS['3bme_type'] = "no2b"
#     ARGS['file2e1max'] = '16 file2e2max=32 file2lmax=16'
#     ARGS['2bme'] = '/global/scratch/exch/me2j/takayuki/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw%d_emax16_e2max32.me2j.gz'%(hw)
#     ARGS['file3e1max'] = '16 file3e2max=32 file3e3max=24'
#     ARGS['3bme'] = '/global/scratch/exch/me3j/takayuki/NO2B_compact/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms16_32_24.stream.bin'%(hw)
#     ARGS['3bme_type'] = "no2b"

#     ARGS['LECs'] = 'N2LO_sat'
##     ARGS['LECs'] = 'n7_N2LO_sat'
##     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
##     ARGS['2bme'] = '/global/scratch/exch/ME_share/TBMEA2n2losat_14_28.%d_TUD.int.gz'%(hw)
##     ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
##     ARGS['3bme'] = '/global/scratch/exch/ME_share/v3trans_J3T3.int_NN3Nnnlosat_nu3_330_161615.%d_form.gz'%(hw)
#     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
#     ARGS['2bme'] = '/global/scratch/exch/ME_share/TBMEA2n2losat_14_28.%d_TUD.int.gz'%(hw)
#     ARGS['file3e1max'] = '18 file3e2max=36 file3e3max=24'
#     ARGS['3bme'] = '/global/scratch/exch/me3j/takayuki/NO2B_compact/NO2B_ThBME_N2LOsat_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin'%(hw)
#     ARGS['3bme_type'] = "no2b"

#     ARGS['LECs'] = 'DNNLOgo'
#     ARGS['file2e1max'] = '18 file2e2max=36 file2lmax=18'
#     ARGS['2bme'] = '/global/scratch/exch/me2j/takayuki/TwBME-HO_NN-only_DNNLOgo_bare_hw%d_emax18_e2max36.me2j.gz'%(hw)
##     ARGS['file3e1max'] = '18 file3e2max=36 file3e3max=22'
##     ARGS['3bme'] = '/global/scratch/exch/me3j/takayuki/NO2B_ThBME_DNNLOgo_3NFJmax15_IS_hw%d_ms18_36_22.stream.bin'%(hw)
#     ARGS['file3e1max'] = '16 file3e2max=32 file3e3max=26'
#     ARGS['3bme'] = '/global/scratch/exch/me3j/takayuki/NO2B_compact/NO2B_ThBME_DNNLOgo_3NFJmax15_IS_hw%d_ms16_32_26.stream.bin'%(hw)
#     ARGS['3bme_type'] = "no2b"

#     ARGS['LECs'] = 'N3LO_3N400'
#     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
#     ARGS['2bme'] = '/global/scratch/exch/ME_share/TBMEA2n3lo-srg2.0_14_28.%d_TUD.int.gz'%(hw)
#     ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
#     ARGS['3bme'] = '/global/scratch/exch/ME_share/v3trans_J3T3.int_N3LO3NF400_-0.2-srg2.0_330_161615.%d_form.gz'%(hw)

# Corrected N4LO LNL potential
#     ARGS['LECs'] = 'N4LO_LNL2'
##     ARGS['LECs'] = 'OccFile_N4LO_LNL2'
#     ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
#     ARGS['2bme'] = '/global/scratch/exch/ME_share/TBMEA2n4lo500-srg2.0_14_28.%d_TUD.int.gz'%(hw)
##     ARGS['file3e1max'] = '18 file3e2max=36 file3e3max=22'
##     ARGS['3bme'] = '/home/holt/projects/deft
##     holt/tmiyagi/MtxElmnt/3BME/NO2B_ThBME_srg2.00_Nmax40_N4LO_EMN500_ChEFT_N2LO_cD-1.80cE-0.31_LNL2_IS_hw16_ms18_36_22.stream.bin'%(hw)
##     ARGS['3bme_type'] = "no2b"
#     ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
#     ARGS['3bme'] = '/global/scratch/exch/ME_share/v3trans_J3T3.int_NNn4lo500_3NlnlcD-1.8cE-0.31_-1.8-srg2.0_220_161616.%d_from24_form.gz'%(hw)

#     ARGS['LECs'] = 'N4LO_srg2.0'
#     ARGS['2bme'] = '/global/scratch/exch/me2j/takayuki/TwBME-HO_NN-only_N4LO_EMN500_srg2.00_hw%d_emax6_e2max12.me2j.gz'%(hw)
#     ARGS['3bme'] = 'none'

#     ARGS['LECs'] = 'srg0625'
#     ARGS['2bme'] = '/global/home/bhu/IMSRG/EOM-IMSRG/TBME_inputs/chi2b_srg0625_eMax08_hwHO0%d.me2j.gz'%(hw)
#     ARGS['2bme'] = '/global/scratch/exch/me2j/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
#     ARGS['3bme'] = '/global/scratch/exch/me3j/chi2b3b400cD-02cE0098_hwconv036_srg0625ho40J_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)

#     ARGS['LECs'] = 'srg0800'
#     ARGS['2bme'] = '/global/scratch/bhu/TwBME-HO_NN-only_N3LO_EM500_srg1.88_hw%d_emax4_e2max8.me2j.gz'%(hw)

     ARGS['LECs'] = 'usdb'
     ARGS['fmt2'] = 'nushellx'
     ARGS['2bme'] = '../input/usdbpn_bhu.int'
##     ARGS['fmt2'] = 'tokyo'
##     ARGS['2bme'] = '../input/usdb.snt'
     ARGS['3bme'] = 'none'
     ARGS['basis'] = 'oscillator'

#     ARGS['LECs'] = 'srg0625'
#     ARGS['2bme'] = '/work/hda21/hda215/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
#     ARGS['3bme'] = '/work/hda21/hda215/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
#     ARGS['2bme'] = '/global/scratch/exch/me2j/chi2b_srg0625_eMax14_lMax10_hwHO0%d.me2j.gz'%(hw)
#     ARGS['3bme'] = '/itch/exch/me3j/new/chi2b3b400cD-02cE0098_hwconv036_srg0625ho40J_eMax14_EMax14_hwHO0%d.me3j.gz'%(hw)

     ARGS['hw'] = '%d'%hw
     #ARGS['hw'] = '%.4f'%hw
     ARGS['A'] = '%d'%A
#     ARGS['valence_space'] = 'GSM'
#     ARGS['valence_space'] = reference
     ARGS['valence_space'] = 'sd-shell'
#     ARGS['valence_space'] = 'psd-shell'
#     ARGS['valence_space'] = '0hw-shell'
#     ARGS['valence_space'] = 'sdfp-shell'
#     ARGS['valence_space'] = 'spsdpf-shell'
#     ARGS['valence_space'] = 'Cr%d'%A
#     ARGS['valence_space'] = 'psd5-He8'
#     ARGS['custom_valence_space'] = 'He8,p0p1,p0p3,n0p1,n0d5,n1s1'
#     ARGS['valence_space'] = 'psd5-n-shell'
#     ARGS['custom_valence_space'] = 'He4,p0p1,p0p3,n0p1,n0p3,n0d5,n1s1'
#     ARGS['valence_space'] = 'psd5-shell'
#     ARGS['custom_valence_space'] = 'He4,p0p1,p0p3,p0d5,p1s1,n0p1,n0p3,n0d5,n1s1'
#     ARGS['valence_space'] = 'sdf7p-O16'
#     ARGS['custom_valence_space'] = 'O16,p0d5,n0d5,p0d3,n0d3,p1s1,n1s1,n0f7,n1p3,n1p1'
#     ARGS['valence_space'] = 'sdf7p-O22'
#     ARGS['custom_valence_space'] = 'O22,p0d5,p0d3,n0d3,p1s1,n1s1,n0f7,n1p3,n1p1'
#     ARGS['valence_space'] = 'sdf7p-O24'
#     ARGS['custom_valence_space'] = 'O24,p0d5,p0d3,p1s1,n0d3,n0f7,n1p3,n1p1'
#     ARGS['valence_space'] = 'sdfpg9-36O'
#     ARGS['custom_valence_space'] = 'O36,p0d5,p0d3,p1s1,n0f5,n1p3,n1p1,n0g9'
#     ARGS['valence_space'] = 'pf5g9-Ni56'
#     ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9'
#     ARGS['valence_space'] = 'sdg7h11-Sn100'
#     ARGS['custom_valence_space'] = 'Sn100,p0g7,n0g7,p1d5,n1d5,p1d3,n1d3,p2s1,n2s1,p0h11,n0h11'

#     ARGS['Operators'] = ''    # Operators to consistenly transform, separated by commas.
#     ARGS['Operators'] = 'Rp2'
#     ARGS['Operators'] = 'Rp2,Rn2'
#     ARGS['Operators'] = 'E2'
#     ARGS['Operators'] = 'M1,E2,Rp2,Rn2'
#     ARGS['Operators'] = 'Rp2,Rn2,IVD'
#     ARGS['Operators'] = 'Rp2,Rn2,Rm2,ISM'
#     ARGS['Operators'] = 'Rp2,Rn2,Rm2,GamowTeller'
#     ARGS['Operators'] = 'Rp2,Rn2,Rm2,IVD,ISM'
#     ARGS['Operators'] = 'E2,M1,GamowTeller'
#     ARGS['Operators'] = 'M1p,M1n,Sigma_p,Sigma_n'
#     ARGS['Operators'] = 'GamowTeller'

     if ARGS['LECs'] == 'srg0800' : c3='-3.2'; c4='5.4'; cD='0.0'
     #elif ARGS['LECs'] == 'usdb' : c3='-3.2'; c4='5.4'; cD='-1.8'
     elif ARGS['LECs'] == 'usdb' : c3='-3.2'; c4='5.4'; cD='0.0'
     elif ARGS['LECs'] == 'EM1.8_2.0' : c1='-0.81'; c3='-3.2'; c4='5.4'; cD='1.264'
     #elif ARGS['LECs'] == 'EM1.8_2.0' : c1='-0.81'; c3='-3.2'; c4='5.4'; cD='0.0'
     elif ARGS['LECs'] == 'N4LO_LNL2' : c1='-1.1'; c3='-5.54'; c4='4.17'; cD='-1.8'
     #elif ARGS['LECs'] == 'N4LO_LNL2' : c1='-1.1'; c3='-5.54'; c4='4.17'; cD='0.0'
     elif ARGS['LECs'] == 'DNNLOgo' : c1='-0.74'; c3='-0.65'; c4='0.96'; cD='0.081'
     elif ARGS['LECs'] == 'N2LO_sat' : c1='-1.1215'; c3='-3.9250'; c4='3.7657'; cD='0.8168'
     else: print("No default c3 and c4 for interaction " + ARGS['LECs']); break;

     op = ['Long5','ET'] #['Long5', 'ET', 'MT']
     #op = ['Long5'] #['Long5', 'ET', 'MT']
     rho = ['0.120'] #['0.100','0.105','0.110','0.115','0.120']
     #rho = ['0.100','0.105','0.110','0.115','0.120']
     #rho = ['0.100','0.105','0.110','0.115','0.120'] #['0.100','0.105','0.110','0.115','0.120']
     #struct_facts = ['p'] #['00','p','n']
     struct_facts = ['00','p','n'] #['00','p','n']
     op_indexes = ['1'] #['0', '1']
     L = ['1'] #['1', '3', '5']
     #L = ['1', '3', '5']
     #p = ['433.34']
     #p = ['0.0','25','50.0','75.0','125.0','150.0','159.98473','175.0','225.02806','275.00894','314.98507','354.98784'] # for S00
     #p = ['0.00','86.73','151.26','195.55','231.51','262.59','267.32','271.97','274.27','281.04','285.47','289.83','312.72','334.04','354.08'] # for Sp
     #p = ['0.00','106.15','185.12','239.32','283.33','312.72','321.37','327.16','332.85','338.66','343.95','349.37','354.70','360.00','375.00','393.00','433.34','484.00'] # for Sp of F19
     #p = ['0.00','106.15','185.12','239.32','283.33','321.37','327.16','332.85','335.66','343.95','349.37','354.70','382.71','408.81','433.34'] # for Sn
     #p = [0.0,25,50.0,75.0,112.1,125.0,137.3,150.0,159.98473,175.0,177.2,194.1,209.7,225.02806, \
     #     237.7,250.6,262.8,275.00894,285.7,296.5,306.9,314.98507,326.7,336.2,345.4,354.98784]
     if(ARGS['LECs']=='usdb'):
         c3s = ['-3.20','-2.20','-3.40','-2.40','-4.78','-3.78','0.00']
         c4s = ['5.40','4.40','3.40','2.40','3.96','2.96','0.00']
     else:
         c3s = [c3]
         c4s = [c4]
     DMSDEFT=''
     for iop in range(0,len(op)):
         for isf in range(0,len(struct_facts)):
             for ir in range(0,len(rho)):
                 if struct_facts[isf]=='00' and ir >0 : continue
                 for id in range(0,len(op_indexes)):
                     for iL in range(0,len(L)):
                         for ics in range(0,len(c3s)):
                             c3 = c3s[ics]; c4 = c4s[ics]
                             #p = ['0.00']
                             p = ['0.00']
                             #p = []
                             #fobj = open('momenta_list_S'+struct_facts[isf]+'_'+ARGS['reference']+'_hw'+ARGS['hw']+'_A'+ARGS['A']+'.dat','r')
                             #fobj = open('momenta_list_'+ARGS['reference']+'_A'+ARGS['A']+'.dat','r')
                             #for eachline in fobj: p.append(eachline.strip('\n'))
                             #fobj.close()
                             for ip in range(0,len(p)):
                                 DMSDEFT += 'DMSDEFT_'+op[iop]+'_'+rho[ir]+'_'+c3+'_' \
                                 +c4+'_'+cD+'_'+struct_facts[isf]+'_'+L[iL]+'_'+p[ip]+','
     #ARGS['Operators'] = DMSDEFT[:-1] + ',Rp2,Sigma_p,Sigma_n'
     ARGS['Operators'] = 'Sigma_p,Sigma_n'
     #ARGS['Operators'] = 'Rp2,Rn2'
     #ARGS['Operators'] = 'Rp2,E2'
     #ARGS['Operators'] = DMSDEFT[:-1]
     #print(ARGS['Operators'])


    ### Make an estimate of how much time to request. Only used for slurm at the moment.
     time_request = '120:00:00'
     if   e <  5 : time_request = '00:10:00'
     elif e <  8 : time_request = '24:00:00'
     elif e < 10 : time_request = '48:00:00'
     elif e < 12 : time_request = '96:00:00'

#     jobname  = '%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['smax'],ARGS['hw'],ARGS['A'])
     jobname  = '%s_%s_%s_%s_e%s_E%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['hw'],ARGS['A'])

  ### Some optional parameters that we probably want in the output name if we're using them
     if 'lmax3' in ARGS:  jobname  += '_l%d'%(ARGS['lmax3'])
     if 'eta_criterion' in ARGS: jobname += '_eta%s'%(ARGS['eta_criterion'])
     if 'core_generator' in ARGS: jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS: jobname += '_' + ARGS['BetaCM']
     if 'denominator_delta' in ARGS: jobname += '_dE' + ARGS['denominator_delta']
     if 'emax_imsrg' in ARGS: jobname += '_eimsrg' + ARGS['emax_imsrg']
     if 'e2max_imsrg' in ARGS: jobname += '_e2imsrg' + ARGS['e2max_imsrg']
     if 'basis' in ARGS: jobname += '_' + ARGS['basis']
     ARGS['flowfile'] = 'output/BCH_' + jobname + '.dat'
     ARGS['intfile']  = 'output/' + jobname
#     ARGS['intfile']  = '/global/scratch/bhu/output/' + jobname

     logname = jobname + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')

     cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

  ### Submit the job if we're running in batch mode, otherwise just run in the current shell
     if batch_mode==True:
       sfile = open(jobname+'.batch','w')
       if BATCHSYS == 'PBS':
         sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,mail_address,logname,NTHREADS,cmd))
         sfile.close()
         call(['qsub', jobname+'.batch'])
       elif BATCHSYS == 'SLURM':
         sfile.write(FILECONTENT%(NTHREADS,jobname,time_request,mail_address,NTHREADS,NTHREADS,cmd))
         sfile.close()
         call(['sbatch', jobname+'.batch'])
       remove(jobname+'.batch') # delete the file
       sleep(0.1)
     else:
       call(cmd.split())  # Run in the terminal, rather than submitting

