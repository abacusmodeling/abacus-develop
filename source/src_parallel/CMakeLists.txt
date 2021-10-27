add_library(
    parallel
    OBJECT
    ft.cpp
    parallel_common.cpp
    parallel_global.cpp
    parallel_grid.cpp
    parallel_kpoints.cpp
    parallel_orbitals.cpp
    parallel_pw.cpp
    parallel_reduce.cpp
    subgrid_oper.cpp
)

if(ENABLE_DEEPKS)
	target_sources(parallel PRIVATE parallel_deepks.cpp)
endif()

