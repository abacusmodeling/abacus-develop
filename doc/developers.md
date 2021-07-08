# For developers

- [Build and install ABACUS with CMake](#build-and-install-abacus-with-cmake)
- [Raising issues on GitHub](#raising-issues-on-github)
- [Modularization and module tests](#modularization-and-module-tests)
- [Contributing to ABACUS](#contributing-to-abacus)
    - [Making pull requests](#making-pull-requests)
    - [Providing tests](#providing-tests)
    - [Upating documentation](#updating-documentation)

[back to main page](../README.md)

## Build and install ABACUS with CMake

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

[back to top](#for-developers)

## Raising issues on GitHub
Raise issues on GitHub is a convernient way to notify the develper team about bugs and feature requests of the ABACUS code. We provide a few templates for issues.

[back to top](#for-developers)

## Modularization and module tests
The ABACUS code is reconstructed to form several self-contained modules. A description of modules can be found in the [installation guide](install.md#structure-of-source-code). We also provide module tests for the modules.

[back to top](#for-developers)

## Contributing to ABACUS

### Making pull requests

### Providing tests

### Updating documentation

[back to top](#for-developers)