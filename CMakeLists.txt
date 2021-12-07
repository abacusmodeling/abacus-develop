########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.15)
########################################

project(ABACUS
    VERSION 2.2.0
    DESCRIPTION "ABACUS is an electronic structure package based on DFT."
    HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
    # LANGUAGES CXX
)

option(ENABLE_DEEPKS "Enable DeePKS functionality" OFF)
option(ENABLE_LIBXC "Enable LibXC functionality" OFF)
option(USE_CUDA "Enable support to CUDA." OFF)
option(USE_ROCM "Enable support to ROCm." ON)
option(USE_OPENMP " Enable OpenMP in abacus." ON)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(BUILD_TESTING "Build ABACUS unit tests" OFF)
option(GENERATE_TEST_REPORTS "Enable test report generation" OFF)

set(ABACUS_BIN_NAME abacus)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(ABACUS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(ABACUS_BIN_PATH ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME})
include_directories(${ABACUS_SOURCE_DIR})
add_executable(${ABACUS_BIN_NAME} source/main.cpp)

find_package(Cereal REQUIRED)
include_directories(${Cereal_INCLUDE_DIR})
add_compile_definitions(USE_CEREAL_SERIALIZATION)

find_package(ELPA REQUIRED)
include_directories(${ELPA_INCLUDE_DIR})
target_link_libraries(${ABACUS_BIN_NAME} ELPA::ELPA)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(${ABACUS_BIN_NAME} MPI::MPI_CXX)
add_compile_definitions(__MPI)

find_package(Threads REQUIRED)
target_link_libraries(${ABACUS_BIN_NAME} Threads::Threads)

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(${ABACUS_BIN_NAME} OpenMP::OpenMP_CXX)
  # add_compile_options(${OpenMP_CXX_FLAGS})
  add_compile_definitions(__OPENMP)
endif()

set(CMAKE_CXX_STANDARD 11)
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  if(NOT DEFINED USE_CUDA)
    message("CUDA components detected. \nWill build the CUDA version of ABACUS.")
    set(USE_CUDA ON)
  else()
    if(NOT USE_CUDA)
      message(WARNING "CUDA components detected, but USE_CUDA set to OFF. \nNOT building CUDA version of ABACUS.")
    endif()
  endif()
else() # CUDA not found
  if (USE_CUDA)
    message(FATAL_ERROR "USE_CUDA set but no CUDA components found.")
    set(USE_CUDA OFF)
  endif()
endif()
if(USE_CUDA)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_EXTENSIONS ON)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(CUDA)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lcufft
    -lcublas
    -lcudart
    -lcusolver
  )
  set_property(TARGET ${ABACUS_BIN_NAME}
  PROPERTY CUDA_ARCHITECTURES
  60 # P100
  70 # V100
  75 # T4
  )
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  add_compile_definitions(__CUDA)
endif()

if(USE_ROCM)
  set(CMAKE_CXX_STANDARD 11) 
  SET(CMAKE_HIP_STANDARD 11)
  # set(CMAKE_CXX_STANDARD_REQUIRED True) 
  list(APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
  list(APPEND CMAKE_SYSTEM_PREFIX_PATH "/opt/rocm/hip")
  include_directories(
    /opt/rocm/hip/include
    /opt/rocm/hipfft/include
    /opt/rocm/hipblas/include
    /opt/rocm/rocsolver/include
  )
  # find_package(ROCM REQUIRED)
  find_package(hip REQUIRED)
  find_package(hipfft REQUIRED)
  find_package(hipblas REQUIRED)
  find_package(rocsolver REQUIRED)
  target_link_libraries(${ABACUS_BIN_NAME}
    hip::device
    hip::host
    hip::hipfft
    roc::hipblas
    roc::rocsolver
  )
  add_compile_definitions(__ROCM)
  # set (HIP_HIPCC_FLAGS -hc; -fno-gpu-rdc; --amdgpu-target=gfx906)
  # add_compile_definitions(__HIP_PLATFORM_HCC__)
endif()

if(ENABLE_ASAN)
  add_compile_options(
    -O1
    -g
    -fsanitize=address
    -fno-omit-frame-pointer
  )
  add_link_options(
    -g
    -fsanitize=address
  )
else()
  add_compile_options(-O2 -g)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings" )
endif()

if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
endif()
if(MKLROOT)
  find_package(IntelMKL REQUIRED)
  add_definitions(-D__MKL -DMKL_ILP64)
  include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lifcore
    IntelMKL::MKL
  )
else()
  find_package(FFTW3 REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(ScaLAPACK REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIRS})
  target_link_libraries(${ABACUS_BIN_NAME}
    FFTW3::FFTW3
    LAPACK::LAPACK
    ScaLAPACK::ScaLAPACK
  )
endif()

if(ENABLE_DEEPKS)
  set(CMAKE_CXX_STANDARD 14)
  find_package(Torch REQUIRED)
  include_directories(${TORCH_INCLUDE_DIRS})
  target_link_libraries(${ABACUS_BIN_NAME} ${TORCH_LIBRARIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

  find_path(libnpy_SOURCE_DIR
    npy.hpp
    HINTS ${libnpy_INCLUDE_DIR}
  )
  if(NOT libnpy_SOURCE_DIR)
  include(FetchContent)
    FetchContent_Declare(
      libnpy
      GIT_REPOSITORY https://github.com.cnpmjs.org/llohse/libnpy.git
    )
    FetchContent_MakeAvailable(libnpy)
  endif()
  include_directories(${libnpy_SOURCE_DIR})
  add_compile_definitions(__DEEPKS)
endif()

if(DEFINED Libxc_DIR)
  set(ENABLE_LIBXC ON)
endif()
if(ENABLE_LIBXC)
  find_package(Libxc)
  if(${Libxc_FOUND})
    message("Using Libxc.")
    add_compile_definitions(USE_LIBXC)
    target_link_libraries(${ABACUS_BIN_NAME} Libxc::xc)
    include_directories(${Libxc_INCLUDE_DIRS})
  else()
    message(WARNING "Will not use Libxc.")
  endif()
endif()

add_compile_definitions(
  __EXX
  __FFTW3
  __FP
  __SELINV
  __LCAO
  METIS
  EXX_DM=3
  EXX_H_COMM=2
  TEST_EXX_LCAO=0
  TEST_EXX_RADIAL=1
)

add_subdirectory(source)

target_link_libraries(${ABACUS_BIN_NAME}
    base
    cell
    symmetry
    md
    symmetry
    neighbor
    orb
    io
    ions
    lcao
    parallel
    mrrr
    pdiag
    pw
    ri
    driver
    -lm
)
if(CMAKE_COMPILER_IS_GNUCXX)
target_link_libraries(${ABACUS_BIN_NAME}
  -lgfortran
)
endif()

install(PROGRAMS ${ABACUS_BIN_PATH}
    TYPE BIN
    #DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CTest)
enable_testing()
IF (BUILD_TESTING)
  add_subdirectory(tests)
endif()
