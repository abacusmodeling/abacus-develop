list(APPEND objects
    xc_functional.cpp
    xc_functional_vxc.cpp
    xc_gga_pw.cpp
    xc_funct_exch_lda.cpp
    xc_funct_corr_lda.cpp
    xc_funct_exch_gga.cpp
    xc_funct_corr_gga.cpp
    xc_funct_hcth.cpp
    )

add_library(
    xc
    OBJECT
    ${objects}
)
