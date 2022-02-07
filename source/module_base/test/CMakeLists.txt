remove_definitions(-D__MPI)
AddTest(
  TARGET base_tool_quit
  SOURCES tool_quit_test.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_tool_check
  SOURCES tool_check_test.cpp ../tool_check.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_memory
  SOURCES memory_test.cpp ../memory.cpp
)
AddTest(
  TARGET base_matrix3
  LIBS ${math_libs}
  SOURCES matrix3_test.cpp ../matrix3.cpp ../matrix.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_blas_connector
  LIBS ${math_libs}
  SOURCES blas_connector_test.cpp
)
AddTest(
  TARGET base_complexarray
  SOURCES complexarray_test.cpp ../complexarray.cpp
)
AddTest(
  TARGET base_integral
  SOURCES math_integral_test.cpp
)