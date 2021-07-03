########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.18)
########################################

project(ABACUS
    VERSION 2.2.0
    DESCRIPTION "ABACUS is an electronic structure package based on DFT."
    HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
    LANGUAGES CXX
)
set(ABACUS_BIN_NAME ${PROJECT_NAME}-${PROJECT_VERSION})

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-O2 -g)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)

set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${ABACUS_SOURCE_DIR})

# We only use headers
find_package(Boost REQUIRED)
find_package(Cereal REQUIRED)

include_directories(
    ${Cereal_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# We only use libraries while headers are customed in project:
# - "blas_connector.h"
# - "lapack_connect.h"
# - "my_elpa.h"
# - "scalapack_connector.h"
find_package(BLAS REQUIRED)
find_package(ELPA REQUIRED)
find_package(ScaLAPACK REQUIRED)

# Both header and library are needed
find_package(FFTW3 REQUIRED)
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH} ${FFTW3_INCLUDE_DIRS})

# No header needed
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

add_definitions(
    -D__EXX
    -D__FFTW3
    -D__FP
    -D__MPI
    -D__OPENMP
    -D__SELINV
    -D__LCAO
    -DMETIS
    -DEXX_DM=3
    -DEXX_H_COMM=2
    -DTEST_EXX_LCAO=0
    -DTEST_EXX_RADIAL=1
    -DUSE_CEREAL_SERIALIZATION
)

add_subdirectory(source)

add_executable(${ABACUS_BIN_NAME}
    source/main.cpp
)

target_link_libraries(${ABACUS_BIN_NAME}
    base
    cell
    neighbor
    orb
    global
    io
    ions
    lcao
    parallel
    mrrr
    pdiag
    pw
    ri
    driver
    -lgfortran
    -lm
    ELPA::ELPA
    FFTW3::FFTW3
    ScaLAPACK::ScaLAPACK
    BLAS::BLAS
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    Threads::Threads
)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME}
    TYPE BIN
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME ${PROJECT_NAME}
)
