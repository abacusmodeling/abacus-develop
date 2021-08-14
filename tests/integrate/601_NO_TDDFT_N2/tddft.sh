#!/bin/bash
#SBATCH -J tddft
#SBATCH -p cn_nl
#SBATCH -N 1
#SBATCH -o tddft.out
#SBATCH -e tddft.err
#SBATCH --no-requeue
#SBATCH -A mhchen_g1
#SBATCH --qos=mhchencnnl
#SBATCH -n 1
hosts=`scontrol show hostname $SLURM_JOB_NODELIST` ;echo $hosts
source /appsnew/source/intel2018.sh
mpirun -n 1 -hosts=${hosts}  /home/mhchen_pkuhpc/mhchen_coe/lustre2/4_liyuanbo/soft/abacus-develop/bin/ABACUS.mpi
