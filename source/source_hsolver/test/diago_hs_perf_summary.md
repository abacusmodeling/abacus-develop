# MPI Performance Summary for PR 7401

## Test setup
- Branch: `pr-7401`
- Binary: `MODULE_HSOLVER_diago_hs_parallel`
- Benchmark driver: `source/source_hsolver/test/diago_hs_perf.sh`
- Problem size: `ndim=240`, `nb=32`, `nbands=120`
- Runs: `case_numb=3`, `loop_numb=3`
- MPI ranks: `1`, `2`, `4`, `8`

## Results
- 1 process:
  - Average Lapack time: `74.33 ms`
  - Average Scalapack time: `85.00 ms`
  - Speedup: `0.87`
- 2 processes:
  - Average Lapack time: `75.67 ms`
  - Average Scalapack time: `82.67 ms`
  - Speedup: `0.92`
- 4 processes:
  - Average Lapack time: `77.67 ms`
  - Average Scalapack time: `153.33 ms`
  - Speedup: `0.51`
- 8 processes:
  - Average Lapack time: `75.67 ms`
  - Average Scalapack time: `81.00 ms`
  - Speedup: `0.93`

## Notes
- The benchmark driver and MPI execution path are working as expected.
- For this test size, distributed `Scalapack` does not beat local `Lapack`.
- The current code is ready for larger-scale MPI testing where speedup may become visible.
