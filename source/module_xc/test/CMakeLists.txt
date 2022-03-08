remove_definitions(-DUSE_LIBXC)

AddTest(
  TARGET XCTest_PBE
  SOURCES test_xc_pbe.cpp ../xc_functional.cpp
)
