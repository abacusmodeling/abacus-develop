add_library(
    plumed_interface
    OBJECT
    plumed_interface.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(plumed_interface)
endif()
