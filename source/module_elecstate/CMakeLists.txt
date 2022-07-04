add_library(
    elecstate
    OBJECT
    elecstate.cpp
    elecstate_pw.cpp
    elecstate_pw_sdft.cpp
    elecstate_lcao.cpp
    dm2d_to_grid.cpp
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
