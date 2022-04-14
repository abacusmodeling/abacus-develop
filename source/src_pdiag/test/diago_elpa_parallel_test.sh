#!/bin/bash

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 2 3 4 6;do
    if [[ $i -gt $np ]];then
        break
    fi
    echo "TEST ELPA in parallel, nprocs=$i"
    mpirun -np $i ./hsolver_diago_elpa
done


