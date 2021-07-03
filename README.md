<p style="text-align:center">
    <img src="https://raw.githubusercontent.com/deepmodeling/abacus-develop/develop/ABACUS.develop/documents/abacus-logo.jpg"/>
</p>

<p style="text-align:center">
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml/badge.svg">
    </a>
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml/badge.svg">
    </a>
</p>

WELCOME TO THE "ABACUS" PROGRAM!

THE PROJECT STARTS FROM https://github.com/abacusmodeling/abacus-develop,
WHERE MORE INFORMATION CAN BE FOUND.

ABACUS IS AN ELECTRONIC STRUCTURE PACKAGE BASED ON DENSITY FUNCTIONAL THEORY.
ABACUS ADOPTS EITHER PLANE WAVE BASIS OR NUMERICAL ATOMIC ORBITALS

---

ABACUS provides the following features and functionalities:

1. Ground-state total energy calculations using Kohn-Sham (KS) density functional theory
(DFT) with local-density, generalized gradient approximations (LDA/GGAs), and hybrid functionals
(PBE0 and HSE06, only for LCAO).
2. Brillouin zone sampling using the Monkhorst-Pack special k-points.
3. Geometry relaxations with Conjugated Gradient (CG) and BFGS methods.
4. Semi-empirical van der Waals energy correction using the Grimme DFT-D2/D3 scheme.
5. NVT molecular dynamics simulation.
6. Stress calculations and cell relaxations.
7. Electric polarization calculation using Berry Phase theory.
8. Interface to the Wannier90 package.
9. Real-time time dependent density functional theory (TDDFT).
10. Electrostatic potential.
11. Mulliken charge analysis.
12. Projected density of states (PDOS).

# Build and install ABACUS from CMake

Check the cmake version on your machine
```bash
cmake --version
```
ABACUS requires the minimum cmake version `3.18`.

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`.
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH}
```
You can provide root path of each dependent package if the package cannot be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDEDIR`, `BOOST_INCLUDEDIR`, `MPI_CXX_COMPILER` and `MKL_DIR`. are currently available to specify.
For example
```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3
```

If the cmake has executed successfully, then
```bash
cmake --build build
cmake --install build
```
If no install prefix is specified, the binary will be installed to `/usr/local/bin/ABACUS` by default.
