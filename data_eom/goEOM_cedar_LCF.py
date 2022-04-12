#!/usr/bin/env python

##########################################################################
##  goUniversal.py
##
##  A python script to run or submit jobs for the common use cases
##  of the EOM code. We check whether there is a pbs or slurm
##  scheduler, assign the relevant input parameters, set names
##  for the output files, and run or submit. 
##  Copy from Ragnar Stroberg
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
NTHREADS=48
#exe = '%s/home/bhu/EOM-IMSRG/src/run_EOM'%(environ['PWD'])
#exe = '/home/bhu/EOM-IMSRG/src/run_EOM'
exe = '/home/b62wong/projects/def-holt/b62wong/NeutronStarStuff/EOM-IMSRG/src/run_EOM'

### Flag to swith between submitting to the scheduler or running in the current shell
#batch_mode=False
batch_mode=True
if 'terminal' in argv[1:]: batch_mode=False

### Don't forget to change this. I don't want emails about your calculations...
mail_address = 'bwong1@triumf.ca'

### This comes in handy if you want to loop over Z
ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}

if BATCHSYS == 'PBS':
  FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q batchmpi
#PBS -d %s
#PBS -l walltime=288:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=60gb
#PBS -m ae
#PBS -M %s
#PBS -j oe
#PBS -o log_eom/%s.o
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/scyld/openmpi/1.10/gnu/lib:/opt/openblas/0.2.14/lib:/opt/boost/1.58.0/lib:/opt/gsl/2.1/lib:/opt/intel/mkl/10.2.2.025/lib/em64t:/opt/intel/mkl/10.2.2.025/lib/64:/opt/hdf5/1.8.15-gcc5.2.0/lib:/opt/gcc/7.5.0/lib64:/opt/gcc/7.5.0/lib:/opt/intel/Compiler/11.1/064/lib/intel64:/opt/anaconda/anaconda2-2.4.0/lib/
%s
  """

elif BATCHSYS == 'SLURM':
  FILECONTENT = """#!/bin/bash
#SBATCH --account=rrg-holt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=log_eom/%s.%%j
#SBATCH --time=%s
#SBATCH --mem=187G
#SBATCH --job-name=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""

### Make a directory for the log files, if it doesn't already exist
if not path.exists('log_eom'): mkdir('log_eom')

#input_path = '/home/bhu/scratch/output'
input_path = '/home/b62wong/scratch/output'
ARGS['method'] = 'magnus'
ARGS['smax'] = '500'

#As=[100,114]
#Zs=[50,50]
As=[208]
Zs=[82]
#As=[40,56,78,90]
#Zs=[20,28,28,40]
#As=[40,56]
#Zs=[20,28]
#As=[78]
#Zs=[28]

#for A in range(48,49):
#for A in [90]:
for i in range(0,len(As)):
 Z=Zs[i]
 A=As[i]
 N=A-Z
 for reference in ['%s%d'%(ELEM[Z],A)]:
  ARGS['reference'] = reference
  ARGS['e3max'] = '24'
  for e in [14]:
   for hw in [12]:
#   for hw in [10,12,14,16]:
     ARGS['A'] = '%d'%A
     ARGS['Z'] = '%d'%Z
     ARGS['N'] = '%d'%N    
     ARGS['hw'] = '%d'%hw
     ARGS['emax'] = '%d'%e
     
     ARGS['valence_space'] = reference
#     ARGS['LECs'] = 'N2LO_sat'
#     ARGS['LECs'] = 'DNNLOgo'
#     ARGS['LECs'] = 'DNLOgo450'
#     ARGS['LECs'] = 'DNNLOgo450'
     ARGS['LECs'] = 'EM1.8_2.0'
#     ARGS['LECs'] = 'EM2.0_2.0'
#     ARGS['LECs'] = 'EM2.2_2.0'
#     ARGS['LECs'] = 'PWA2.0_2.0'
#     ARGS['LECs'] = 'EM2.0_2.5'
     ARGS['Operators'] = 'ISQ' #(used to be IVD)

#     ARGS['core_generator'] = 'atan'
#     ARGS['valence_generator'] = 'shell-model-atan'
#     ARGS['denominator_delta_orbit'] = 'all'
#     ARGS['denominator_delta'] = '10'
     ARGS['BetaCM'] = '5'
#     ARGS['hwBetaCM'] = '12'

### Make an estimate of how much time to request. Only used for slurm at the moment.RGS['core_generator'] = 'atan'a
     #time_request = '168:00:00'
     time_request = '71:59:00'

     jobname  = '%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['smax'],ARGS['hw'],ARGS['A'])
     #jobname  = '%s_%s_%s_e%s_E%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['hw'],ARGS['A'])
     logname = jobname + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')

### Some optional parameters that we probably want in the output name if we're using them
     if 'lmax3' in ARGS:  jobname  += '_l%d'%(ARGS['lmax3'])
     if 'eta_criterion' in ARGS: jobname += '_eta%s'%(ARGS['eta_criterion'])
     if 'core_generator' in ARGS: jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS: jobname += '_' + ARGS['BetaCM']

     emax2 = 10
     emax3 = 20
     
     jobname += '_eimsrg%d_e2imsrg%d_NAT'%(emax2, emax3)
     
     print(jobname)

     initfile = jobname+'_LCF_EOM.ini'
     fx = open('./output/'+initfile,'w')

     fx.write('##########################\n')
     fx.write('#### IMSRG INPUT FILE ####\n')
     fx.write('#### KEEP THIS FORMAT ####\n')
     fx.write('##########################\n')
     fx.write('##########################\n')
     fx.write('# ENTER OUTPUT FILE PREFIX\n')
#     fx.write(jobname + '_TDA_' + ARGS['Operators'] + '\n') 
     fx.write(jobname + '_LCF500_' + ARGS['Operators'] + '\n') 
#     fx.write(jobname + ARGS['Operators'] + '\n') 
     fx.write('# ENTER HAMILTONIAN FILE NAME\n')
     fx.write(input_path + '/' + jobname + '.snt\n')
     fx.write('# ENTER eMax, lMax\n')
     fx.write(ARGS['emax']+','+ARGS['emax']+'\n')
     fx.write('# ENTER 3B INTERACTION FILE NAME\n')
     fx.write('none\n') 
     fx.write('# ENTER E3Max (enter "0" for no three body force)\n')
     fx.write('0\n') 
     fx.write('# ENTER IMSRG EVOLUTED OPERATOR INPUT FILE NAME\n')
     fx.write(input_path + '/' + jobname + ARGS['Operators'] + '.op\n')
     fx.write('# ENTER HAMILTONIAN TYPE\n')
     fx.write('# 1: T-V - Tcm -Vcm  2: harmonic trap T+U+V  3. T+V\n')
     fx.write('1\n') 
     fx.write('# ENTER HO SPACING hw\n')
     fx.write(ARGS['hw'] +'\n')
     fx.write('# ENTER NUMBER OF PROTONS\n')
     fx.write(ARGS['Z']+'\n')
     fx.write('# ENTER NUMBER OF NEUTRONS\n')
     fx.write(ARGS['N']+'\n')
     fx.write('# ENTER 1 for HF basis\n')
     fx.write('# OR 2 for HO basis\n')
     fx.write('1\n')
     fx.write('# ENTER 1 for magnus method\n')
     fx.write('# or 2 for traditional ode (y/n): "quads", "trips"\n') 
     fx.write('1,\'n\',\'n\'\n') 
     fx.write('# 0: gs only, 1: EOM, 2: TDA\n')
     fx.write('# ENTER 0 for ground state only\n') 
     fx.write('1\n')
     fx.write('# ENTER 1 TO CALCULATE Hcm, 0 otherwise\n')
     fx.write('0\n')
     fx.write('# ENTER 1 TO CALCULATE Rrms, 0 otherwise\n') 
     fx.write('0\n') 
     fx.write('# ENTER other observables or "none"\n') 
     fx.write('none\n')
     fx.write('# EOM file (standard.eom)\n')
     fx.write(ARGS['Operators']+'_LCF.eom\n')            
     fx.write('# Lawson beta value, c.m. hw\n') 
     fx.write('0,'+ARGS['hw']+'\n') 
     fx.write('# .true. for checkpointing (only applicable to magnus)\n' )
     fx.write('.false.\n')
     fx.write('# write normal ordered bare, read normal ordered bare\n')
     fx.write('.false. , .false.\n')                                     
     fx.write('#  write normal ordered decoupled, read normal ordered decoupled\n')
     fx.write('.false. , .false.\n')
     fx.write('#  write omega, read omega\n')
     fx.write('.false. , .false.\n')
     fx.write('#  write human readable\n')
     fx.write('.false.\n')
     fx.write('########################################################\n')
     fx.write('# NOTES \n')
     fx.write('#\n')
     fx.write("# 1. THIS FILE'S NAME SHOULD BE SUPPLIED AS THE ONLY\n") 
     fx.write('# COMMAND ARGUMENT WITH THE EXECUTABLE "./run_IMSRG"\n')
     fx.write('#\n')
     fx.write('# 2. THIS FILE IS READ BY SUBROUTINE "read_main_input_file"\n')
     fx.write('# which is found in "basic_IMSRG.f90" \n')
     fx.write('# \n')
     fx.write('# 3. FILENAMES ARE READ TO COMMON BLOCK "files" \n')
     fx.write('########################################################\n')

     fx.close()

     cmd = '%s %s'%(exe,' '.join(['%s'%(initfile)]) )

     jobnamep = jobname  + '_EOM'
### Submit the job if we're running in batch mode, otherwise just run in the current shell
     if batch_mode==True:
      sfile = open(jobnamep+'.batch','w')
      if BATCHSYS == 'PBS':
       sfile.write(FILECONTENT%(jobnamep,environ['PWD'],NTHREADS,mail_address,logname,NTHREADS,cmd))
       sfile.close()
       call(['qsub', jobnamep+'.batch'])
      elif BATCHSYS == 'SLURM':
       sfile.write(FILECONTENT%(NTHREADS,jobnamep,time_request,jobnamep,mail_address,NTHREADS,NTHREADS,cmd))
       sfile.close()
       call(['sbatch', jobnamep+'.batch'])
      remove(jobnamep+'.batch') # delete the file
      sleep(0.1)
     else:
      call(cmd.split())  # Run in the terminal, rather than submitting

