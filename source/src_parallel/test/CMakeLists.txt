AddTest(
  TARGET ParaCommon
  LIBS MPI::MPI_CXX 
  SOURCES parallel_common_test.cpp ../../module_base/global_variable.cpp ../parallel_common.cpp
)

install(FILES parallel_common_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME parallel_common_test
      COMMAND ${BASH} parallel_common_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
