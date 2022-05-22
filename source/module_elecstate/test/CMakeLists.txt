remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
  TARGET EState_updaterhok_pw
  LIBS ${math_libs} base 
  SOURCES updaterhok_pw_test.cpp ../../src_pw/charge.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/ft.cpp ../../src_pw/pw_basis_old.cpp ../../src_pw/pw_complement.cpp ../../src_pw/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_pw/occupy.cpp  ../../module_elecstate/elecstate_pw.cpp ../../module_elecstate/elecstate.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
