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

option(ENABLE_DEEPKS "Enable DeePKS functionality" OFF)
option(ENABLE_LIBXC "Enable LibXC functionality" OFF)
option(USE_CUDA "Enable support to CUDA." ON)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(BUILD_TESTING "Build ABACUS unit tests" OFF)
option(GENERATE_TEST_REPORTS "Enable test report generation" OFF)

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
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  add_compile_definitions(__CUDA)
endif()

if(ENABLE_DEEPKS)
  set(CMAKE_CXX_STANDARD 14)
else()
  set(CMAKE_CXX_STANDARD 11)
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

set(ABACUS_BIN_NAME abacus)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(ABACUS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(ABACUS_BIN_PATH ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME})

include_directories(${ABACUS_SOURCE_DIR})

find_package(Cereal REQUIRED)
find_package(ELPA REQUIRED)
find_package(MPI REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

if(DEFINED MKL_DIR)
  set(MKL_DIR_CACHE "${MKL_DIR}")
endif()

include(FetchContent)

include_directories(
  ${Cereal_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH}
)

if(ENABLE_DEEPKS)
  find_package(Torch REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
  FetchContent_Declare(
    libnpy
    GIT_REPOSITORY https://github.com.cnpmjs.org/llohse/libnpy.git
  )
  FetchContent_MakeAvailable(libnpy)
  include_directories(
    ${libnpy_SOURCE_DIR}
    ${TORCH_INCLUDE_DIRS}
  )
  add_definitions(-D__DEEPKS)
endif()

if(DEFINED Libxc_DIR)
  set(ENABLE_LIBXC ON)
endif()
if(ENABLE_LIBXC)
  find_package(Libxc)
  if(${Libxc_FOUND})
    message("Using Libxc.")
    add_compile_definitions(USE_LIBXC)
    # applying -DUSING_Libxc definition
    # target_compile_definitions(${ABACUS_BIN_NAME} PRIVATE $<TARGET_PROPERTY:Libxc::xc,INTERFACE_COMPILE_DEFINITIONS>)
    include_directories(${Libxc_INCLUDE_DIRS})
  else()
    message(WARNING "Will not use Libxc.")
  endif()
endif()

if(DEFINED MKL_DIR_CACHE)
  find_package(IntelMKL REQUIRED)
  add_definitions(-D__MKL -DMKL_ILP64)
  include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
else()
  find_package(FFTW3 REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(ScaLAPACK REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIRS})
endif()

add_compile_definitions(
  __EXX
  __FFTW3
  __FP
  __MPI
  __OPENMP
  __SELINV
  __LCAO
  METIS
  EXX_DM=3
  EXX_H_COMM=2
  TEST_EXX_LCAO=0
  TEST_EXX_RADIAL=1
  USE_CEREAL_SERIALIZATION
)

add_subdirectory(source)

add_executable(${ABACUS_BIN_NAME}
    source/main.cpp
)

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
    ELPA::ELPA
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    Threads::Threads
)

if(DEFINED MKL_DIR_CACHE)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lifcore
    IntelMKL::MKL
    OpenMP::OpenMP_CXX
  )
else()
  target_link_libraries(${ABACUS_BIN_NAME}
    FFTW3::FFTW3
    ScaLAPACK::ScaLAPACK
    LAPACK::LAPACK
    OpenMP::OpenMP_CXX
  )
endif()

# if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
target_link_libraries(${ABACUS_BIN_NAME}
  -lgfortran
)
# endif()

if(ENABLE_DEEPKS)
  target_link_libraries(${ABACUS_BIN_NAME}
    ${TORCH_LIBRARIES}
  )
endif()

if(USE_CUDA)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lcufft
    -lcublas
    -lcudart
  )
  set_property(TARGET ${ABACUS_BIN_NAME}
  PROPERTY CUDA_ARCHITECTURES
  60 # P100
  70 # V100
  75 # T4
  )
endif()

if(${Libxc_FOUND})
  target_link_libraries(${ABACUS_BIN_NAME} Libxc::xc)
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
