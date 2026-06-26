#!/bin/bash
set -e

cd "$(dirname "$0")"

BINARY=./MODULE_HSOLVER_diago_hs_parallel
if [ ! -x "$BINARY" ]; then
    BINARY="../../../build-pr-7401/source/source_hsolver/test/MODULE_HSOLVER_diago_hs_parallel"
fi
if [ ! -x "$BINARY" ]; then
    echo "ERROR: MODULE_HSOLVER_diago_hs_parallel not found in $(pwd) or build-pr-7401 path"
    exit 1
fi

MPI_RUN=mpirun
if mpirun --help 2>&1 | grep -q -- '--allow-run-as-root'; then
    MPI_RUN='mpirun --allow-run-as-root'
fi

for n in 1 2 4 8; do
    echo "============================================================"
    echo "MPI benchmark: $n processes"
    echo "============================================================"
    OMP_NUM_THREADS=1 $MPI_RUN -np "$n" "$BINARY" --perf | tee perf_${n}.log
    echo ""
done
