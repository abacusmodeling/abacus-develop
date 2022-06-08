add_library(
    planewave
    OBJECT
    fft.cpp
    pw_basis.cpp 
    pw_basis_k.cpp 
    pw_distributeg.cpp 
    pw_distributeg_method1.cpp 
    pw_distributeg_method2.cpp
    pw_distributer.cpp 
    pw_init.cpp 
    pw_transform.cpp 
    pw_transform_k.cpp 
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
