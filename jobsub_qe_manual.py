#!/usr/bin/env python

import os, string
import numpy as np

def convert2qe():
    """ Convert pw.inp to quantum espresso form (to submit without ase).
    Changes 'ase3' to 'bfgs' .
    """
    if os.path.isfile(os.path.join(os.getcwd(),'pw.inp')):
        os.system("sed 's/ase3/bfgs/g' pw.inp > tmp.inp ; mv tmp.inp pw.inp ;")
    return

def convert2ase():
    """ Convert pw.inp to ase form (to submit with ase).
    Changes 'bfgs' to 'ase3' .
    """
    if os.path.isfile(os.path.join(os.getcwd(),'pw.inp')):
        os.system("sed 's/bfgs/ase3/g' pw.inp > tmp.inp ; mv tmp.inp pw.inp ;")
    return
#########################################################
#########################################################

### SUBMISSION PARAMETERS ########
queue   = 'suncat-test'
jobfile = 'runqe.sh'
logfile = 'test.log'

cwd = os.getcwd()
repath = '/nfs/slac/g/suncatfs/'
cwd = os.path.join(repath,cwd.split('/a/suncatfs1/u1/')[-1])

jobpath = os.path.join(cwd,jobfile)
logpath = os.path.join(cwd,logfile)
executablepath = os.path.join(repath,'espresso/bin/pw.x')

##################################
if queue.split('-')[0] == 'suncat2':
    nnodes = 12
if queue.split('-')[0] == 'suncat':
    nnodes = 16
#########################################################
#########################################################
subcmd  = 'bsub -n '+str(nnodes)+' -o '+logpath+' -q '+queue+' '+jobpath
print "submitting with "+subcmd

file = open(jobfile,'w')
file.write("""#!/bin/sh
. /nfs/slac/g/suncatfs/vossj/bashspresso.sh
JOBTMP=/scratch/JOB$LSB_BATCH_JID
perHostMpiExec mkdir $JOBTMP
#add -r to copy if subdirectories are needed to setup job
perHostMpiExec cp * $JOBTMP
cd $JOBTMP
perProcMpiExec -wdir $JOBTMP """+executablepath+""" <$LS_SUBCWD/pw.inp
#don't delete *.save if you need e.g. the charge density
rm -r *wfc* *.save
cd ..
cp -r $JOBTMP $LS_SUBCWD
cd $LS_SUBCWD
rm machinefile uniqmachinefile
rm -r $JOBTMP
""")
file.close()  
os.system('chmod +x '+jobfile)

if os.path.isfile(os.path.join(os.getcwd(),'pw.inp')):
    convert2qe()
    if os.path.isfile(executablepath):
        os.system(subcmd)
    else:
        print "incorrect executable path, executable does not exist!"
else:
    print "intialize pw.inp first (generate from ase with espresso(..., onlycreatepwinp=\'filename\') and calc.initialize(atoms))"
