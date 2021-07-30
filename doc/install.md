# Download and install

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Build and install ABACUS with CMake](#build-and-install-abacus-with-cmake)
  - [Build ABACUS with make](#build-abacus-with-make)
- [Structure of the package](#structure-of-the-package)
  - [Structure of source code](#structure-of-source-code)
  [back to main page](../README.md)

## Installation

### Prerequisites

In order to compile ABACUS, users should make sure that the following prerequisites are
present:

- C++ compiler, supporting C++11. For example, [Intel C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/);
- Fortran compiler;
- MPI compiler. The recommended version are [Intel MPI](https://software.intel.com/enus/mpi-library) or [MPICH](https://www.mpich.org/);
- [Boost C++ library](https://www.boost.org/);
- The ScaLAPACK library. For example, [Intel MKL](https://software.intel.com/en-us/mkl)
or [Netlib ScaLAPACK](http://www.netlib.org/scalapack/);
- The [FFTW library](http://www.fftw.org/). ABACUS now supports both FFTW2 and
FFTW3;
- The [ELPA library](https://elpa.mpcdf.mpg.de/);
- The [CEREAL library](https://uscilab.github.io/cereal/);

[back to top](#download-and-install)

### Build and install ABACUS with CMake

We recommend building ABACUS with `cmake` to avoid dependency issues.
ABACUS requires a minimum `cmake` version of `3.18`. Check the version of `cmake`  on your machine with:

```bash
cmake --version
```

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`. If no install prefix is specified, the binary will be installed to `/usr/local/bin/abacus` by default.

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH}
```

You can provide path of each dependent package if the package cannot be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDEDIR`, `BOOST_INCLUDEDIR`, `MPI_CXX_COMPILER` and `MKL_DIR` are currently available to specify.
For example:

```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3 -DBOOST_INCLUDEDIR=/usr/include/boost
```

After configuring, start build and install by:

```bash
cmake --build build
cmake --install build
```

[back to top](#download-and-install)

### Build ABACUS with make

<!-- Before starting to build the program, note that if you are using Intel MKL library, please set the following environmental variable:

```bash
export MKL_NUM_THREAD=1
``` -->

To compile the ABACUS program using legacy `make`, first edit the file `Makefile.vars` under `source` directory:

```bash
cd source/
vi Makefile.vars
```

Specify the location of the compiler and libraries present in your own machine:

```bash
CPLUSPLUS =
CPLUSPLUS_MPI =
FORTRAN =
LAPACK_DIR =
FFTW_DIR =
BOOST_DIR =
ELPA_DIR =
CEREAL_DIR =
```

For example, below is a case where the Intel C++ compiler, Intel MPI are used, along with Intel MKL library. The file Makefile.vars can be set as
follows:

```bash
CPLUSPLUS = icpc
CPLUSPLUS_MPI = mpiicpc
FORTRAN = ifort
LAPACK_DIR = /opt/intel/.../mkl/lib/intel64/
FFTW_DIR = /opt/fftw-3.3.8/
BOOST_DIR = /opt/boost/1.64.0/
ELPA_DIR = /opt/elpa/2016.05.004/
CEREAL_DIR = /opt/cereal/
```

Another example is where GCC, GFORTRAN, MPICH and ScaLAPACK are used:

```bash
CPLUSPLUS = g++
CPLUSPLUS_MPI = mpicxx
FORTRAN = gfortran
SCALAPACK_DIR = /opt/scalapack/
FFTW3_DIR = /opt/fftw-3.3.8/
BOOST_DIR = /opt/boost/1.64.0/
ELPA_DIR = /opt/elpa/2016.05.004/
CEREAL_DIR = /opt/cereal/
```

For this option, it is further required to set the parameter `LIBS` in `Makefile.system`:

```bash
LIBS = \
  -lgfortran -lm \
  -openmp -lpthread \
  ${SCALAPACK_DIR}/lib/libscalapack.a \
  /opt/lapack/lib/liblapack.a \
  /opt/blas/lib/libblas.a \
  /opt/blacs/lib/libblacs.a \
  ${FFTW_LIB} \
  ${ELPA_LIB} \

```

After modifying the `Makefile.vars` file, execute `make` to build the program.

```bash
make -j
```

After the compilation finishes without error messages (except perhaps for some warnings), an executable program `ABACUS.mpi` will be created in directory `bin/`.

[back to top](#download-and-install)

## Structure of the package

Under the ABACUS directory, there are the following subdirectories:

- cmake/

  which contains relevant files for compiling the code with cmake
- documents/

  which contains a copy of the manual in pdf format
- examples/

  which contains some examples
- source/

  which contains the source code and makefiles
- tests/

  which contains test examples
- tools/

  which currently contains the script for generating the numerical atomic orbitals

[back to top](#download-and-install)

### Structure of source code

The source directory further contains the following folders, where the source files of ABACUS are located:

- module_base
- module_cell
- module_grid
- module_md
- module_neighbor
- module_orbital
- obj
- src_external
- src_global
- src_io
- src_ions
- src_lcao
- src_parallel
- src_pdiag
- src_pw
- src_ri

[back to top](#download-and-install)
