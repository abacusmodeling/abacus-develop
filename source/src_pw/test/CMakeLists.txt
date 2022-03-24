remove_definitions(-D__MPI)
AddTest(
  TARGET hsolver_david
  LIBS ${math_libs} base
  SOURCES diago_david_test.cpp ../diago_david.cpp ../../src_parallel/parallel_reduce.cpp
)

install(FILES data-H DESTINATION ${CMAKE_CURRENT_BINARY_DIR})