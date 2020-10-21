#!/bin/bash
#
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -V
ulimit -s unlimited
#source /public/env/env.sh
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
#EXEC=/public/udata/xiaohui/mesia_extra_dm_md2/source/src_lcao/MESIA-ALPHA.fp_mpi.x
EXEC=/public/udata/zhengdy/Git-version/zdy-version-based-on-20180502/abacus-NewGit/ABACUS.1.0.0/bin/ABACUS.mpi.1.0.0
#EXEC=/public/udata/zhengdy/stress-version/stress-test/abacus_dft/ABACUS.1.0.0/bin/ABACUS.mpi.1.0.0
#mpiexec -n 32 -machinefile $PBS_NODEFILE $EXEC > Log.txt
mpiexec -machinefile $PBS_NODEFILE -np $NP $EXEC > stress.txt 
