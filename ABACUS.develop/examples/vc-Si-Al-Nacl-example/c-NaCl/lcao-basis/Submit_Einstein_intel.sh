#!/bin/bash
#
#PBS -l nodes=1:ppn=16
#PBS -l walltime=960:00:00
#PBS -j oe
#PBS -V
ulimit -s unlimited
#source /public/env/env.sh
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
cat $PBS_NODEFILE|uniq
echo $PBS_NODEFILE
echo $LD_LIBRARY_PATH
EXEC=/public/udata/xiaohui/abacus_git/NewGitRepository/20180626-LcaoStressRelax/compile-test/bin/ABACUS.mpi.1.0.0
#EXEC=/public/udata/xiaohui/abacus_git/NewGitRepository/abacus-NewGit/ABACUS.1.0.0/bin/ABACUS.mpi.1.0.0
#EXEC=/public/udata/xiaohui/abacus_git/NewGitRepository/20180619/compile/bin/ABACUS.mpi.1.0.0
/public/intel2017/impi/2017.1.132/bin64/mpirun -machinefile $PBS_NODEFILE -np 8 $EXEC > Log.txt 
