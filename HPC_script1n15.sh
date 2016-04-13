#!/bin/sh
#
#PBS -N JobNameTEST
#PBS -o job.out
#PBS -e job.err
#PBS -m ae
# e-mail is sent if job Aborts or is Executed
#PBS -l walltime=08:00:00
# job time below 12:00:00 -> short cue
#PBS -l vmem=12gb

module load FSL/5.0.6-ictce-5.5.0
module load R

. $FSLDIR/etc/fslconf/fsl.sh

mkdir /user/scratch/gent/vsc418/vsc41855/Tempresultsn15/
mkdir /user/scratch/gent/vsc418/vsc41855/resultsn15/

cd /user/home/gent/vsc418/vsc41855/Effectsize_run1

Rscript HPC_script_run_n15_Files.R

a=1
while [ $a -lt 2 ]
do
Rscript HPC_script_run_n15.R $a
a=`expr $a + 1`
done

rm -R /user/scratch/gent/vsc418/vsc41855/Tempresultsn15/
