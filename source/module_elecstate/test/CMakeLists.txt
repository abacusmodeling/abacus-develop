remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
add_definitions(-D__NORMAL)

AddTest(
  TARGET EState_updaterhok_pw
  LIBS ${math_libs} planewave_serial base_serial
  SOURCES updaterhok_pw_test.cpp ../../src_pw/charge.cpp ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/ft.cpp ../../src_pw/pw_basis_old.cpp ../../src_pw/pw_complement.cpp 
          ../../src_pw/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_pw/occupy.cpp  
          ../../module_elecstate/elecstate_pw.cpp ../../module_elecstate/elecstate.cpp
)
add_library(
    base_serial
    OBJECT
    ../../module_base/complexarray.cpp
    ../../module_base/complexmatrix.cpp
    ../../module_base/element_basis_index.cpp
    ../../module_base/export.cpp
    ../../module_base/integral.cpp
    ../../module_base/inverse_matrix.cpp
    ../../module_base/global_file.cpp
    ../../module_base/global_function.cpp
    ../../module_base/global_function_ddotreal.cpp
    ../../module_base/global_variable.cpp
    ../../module_base/intarray.cpp
    ../../module_base/math_integral.cpp
    ../../module_base/math_polyint.cpp
    ../../module_base/math_sphbes.cpp
    ../../module_base/math_ylmreal.cpp
	  ../../module_base/math_bspline.cpp
    ../../module_base/math_chebyshev.cpp
    ../../module_base/mathzone.cpp
    ../../module_base/mathzone_add1.cpp
    ../../module_base/matrix.cpp
    ../../module_base/matrix3.cpp
    ../../module_base/memory.cpp
    ../../module_base/mymath1.cpp
    ../../module_base/mymath3.cpp
    ../../module_base/polint.cpp
    ../../module_base/realarray.cpp
    ../../module_base/sph_bessel.cpp
    ../../module_base/sph_bessel_recursive-d1.cpp
    ../../module_base/sph_bessel_recursive-d2.cpp
    ../../module_base/timer.cpp
    ../../module_base/tool_check.cpp
    ../../module_base/tool_quit.cpp
    ../../module_base/tool_title.cpp
    ../../module_base/ylm.cpp
)
add_library(
  planewave_serial
  OBJECT
  ../../module_pw/fft.cpp
  ../../module_pw/pw_basis.cpp 
  ../../module_pw/pw_basis_k.cpp 
  ../../module_pw/pw_distributeg.cpp 
  ../../module_pw/pw_distributeg_method1.cpp 
  ../../module_pw/pw_distributeg_method2.cpp
  ../../module_pw/pw_distributer.cpp 
  ../../module_pw/pw_init.cpp 
  ../../module_pw/pw_operation.cpp 
  ../../module_pw/pw_transform.cpp 
  ../../module_pw/pw_transform_k.cpp 
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
