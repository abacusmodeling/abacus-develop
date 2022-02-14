remove_definitions(-D__MPI)
AddTest(
  TARGET base_timer
  SOURCES timer_test.cpp ../timer.cpp
)
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
  SOURCES math_integral_test.cpp ../math_integral.cpp
)
AddTest(
  TARGET base_complexmatrix
  LIBS ${math_libs}
  SOURCES complexmatrix_test.cpp ../complexmatrix.cpp ../matrix.cpp 
)
AddTest(
  TARGET base_matrix
  LIBS ${math_libs}
  SOURCES matrix_test.cpp ../matrix.cpp
)
AddTest(
  TARGET base_sph_bessel_recursive
  SOURCES sph_bessel_recursive_test.cpp ../sph_bessel_recursive-d1.cpp ../sph_bessel_recursive-d2.cpp
)