###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIR - Where to find cereal headers.
#

find_path(Cereal_INCLUDE_DIR
    cereal/cereal.hpp
    HINTS ${CEREAL_INCLUDE_DIR}
)

if(NOT Cereal_INCLUDE_DIR)
    include(FetchContent)
    FetchContent_Declare(
        cereal
        GIT_REPOSITORY https://github.com.cnpmjs.org/USCiLab/cereal.git
        # URL https://github.com/USCiLab/cereal/archive/refs/tags/v1.3.0.tar.gz
    )
    set(Cereal_INCLUDE_DIR ${cereal_SOURCE_DIR})
endif()
# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cereal DEFAULT_MSG Cereal_INCLUDE_DIR)

# Copy the results to the output variables and target.
mark_as_advanced(Cereal_INCLUDE_DIR)
