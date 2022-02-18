list(APPEND objects
    H_XC_pw.cpp
    potential_libxc_meta.cpp
    xc_functional.cpp
    xc_gga_pw.cpp
    )

add_library(
    xc
    OBJECT
    ${objects}
)
