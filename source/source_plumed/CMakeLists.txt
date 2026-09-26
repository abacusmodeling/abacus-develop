add_library(
    plumed_interface
    OBJECT
    plumed_interface.cpp
)

# The PLUMED interface is only compiled in when the build is configured with
# -DPLUMED_ROOT=<prefix of a PLUMED installation>
if(DEFINED PLUMED_ROOT)
  target_compile_definitions(plumed_interface PRIVATE USE_PLUMED)
endif()

if(ENABLE_COVERAGE)
  add_coverage(plumed_interface)
endif()
