remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)

AddTest(
  TARGET EState_updaterhok_pw
  LIBS ${math_libs} base 
  SOURCES updaterhok_pw_test.cpp ../../src_pw/charge.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/ft.cpp ../../src_pw/pw_basis_old.cpp ../../src_pw/pw_complement.cpp ../../src_pw/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_pw/occupy.cpp  ../../module_elecstate/elecstate_pw.cpp ../../module_elecstate/elecstate.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET EState_psiToRho_lcao
  LIBS ${math_libs} ELPA::ELPA base orb cell neighbor
  SOURCES elecstate_lcao_test.cpp ../elecstate_lcao.cpp ../dm2d_to_grid.cpp
            ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp 
            ../../src_lcao/local_orbital_charge.cpp ../../src_lcao/DM_gamma.cpp ../../src_lcao/DM_k.cpp 
            ../../src_lcao/local_orbital_wfc.cpp 
            ../../module_gint/grid_technique.cpp ../../module_gint/grid_meshball.cpp ../../module_gint/grid_bigcell.cpp 
            ../../module_gint/grid_meshcell.cpp ../../module_gint/grid_meshk.cpp
            ../../module_gint/gint_tools.cpp
            ../../src_lcao/dm_2d.cpp
            ../../src_lcao/record_adj.cpp  ../../src_lcao/LCAO_nnr.cpp
            ../../module_gint/gint.cpp ../../module_gint/gint_gamma.cpp
            ../../module_gint/gint_fvl.cpp ../../module_gint/gint_rho.cpp
            ../../module_gint/gint_vl.cpp
            ../../src_pw/charge.cpp 
            ../../src_pdiag/pdiag_common.cpp
            ../../src_io/output.cpp ../../src_pw/soc.cpp
)
target_compile_definitions(EState_psiToRho_lcao PRIVATE __MPI)
install(FILES elecstate_lcao_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME EState_psiToRho_lcao_parallel
      COMMAND ${BASH} elecstate_lcao_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
