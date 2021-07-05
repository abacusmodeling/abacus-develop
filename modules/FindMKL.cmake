# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIRS - where to find mkl.h, etc.
#  MKL_LIBRARIES    - List of libraries when using mkl.
#  MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl_dfti.h HINTS ${MKL_DIR}/include)

find_library(MKL_INTEL NAMES mkl_intel_lp64 HINTS ${MKL_DIR}/lib ${MKL_DIR}/lib/intel64)

find_library(MKL_INTEL_THREAD NAMES mkl_intel_thread HINTS ${MKL_DIR}/lib ${MKL_DIR}/lib/intel64)

find_library(MKL_CORE NAMES mkl_core HINTS ${MKL_DIR}/lib ${MKL_DIR}/lib/intel64)

find_library(MKL_SCALAPACK NAMES mkl_scalapack_lp64 HINTS ${MKL_DIR}/lib ${MKL_DIR}/lib/intel64)

find_library(MKL_BLACS_INTELMPI NAMES mkl_blacs_intelmpi_lp64 HINTS ${MKL_DIR}/lib ${MKL_DIR}/lib/intel64)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI MKL_INCLUDE_DIR)

if(MKL_FOUND)
  set(MKL_LIBRARIES ${MKL_INTEL} ${MKL_INTEL_THREAD} ${MKL_CORE} ${MKL_SCALAPACK} ${MKL_BLACS_INTELMPI})
  set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

  if(NOT TARGET MKL::INTEL)
    add_library(MKL::INTEL UNKNOWN IMPORTED)
    set_target_properties(MKL::INTEL PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET MKL::INTEL_THREAD)
    add_library(MKL::INTEL_THREAD UNKNOWN IMPORTED)
    set_target_properties(MKL::INTEL_THREAD PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL_THREAD}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET MKL::CORE)
    add_library(MKL::CORE UNKNOWN IMPORTED)
    set_target_properties(MKL::CORE PROPERTIES
      IMPORTED_LOCATION "${MKL_CORE}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET MKL::SCALAPACK)
    add_library(MKL::SCALAPACK UNKNOWN IMPORTED)
    set_target_properties(MKL::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${MKL_SCALAPACK}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET MKL::BLACS_INTELMPI)
    add_library(MKL::BLACS_INTELMPI UNKNOWN IMPORTED)
    set_target_properties(MKL::BLACS_INTELMPI PROPERTIES
      IMPORTED_LOCATION "${MKL_BLACS_INTELMPI}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  add_library(MKL::MKL INTERFACE IMPORTED)
  set_property(TARGET MKL::MKL PROPERTY
  INTERFACE_LINK_LIBRARIES
  "-Wl,--start-group"
  MKL::INTEL MKL::INTEL_THREAD MKL::CORE MKL::SCALAPACK MKL::BLACS_INTELMPI
  "-Wl,--end-group"
  )
endif()

mark_as_advanced(MKL_INCLUDE_DIR MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI)
