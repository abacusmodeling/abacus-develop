if(ENABLE_DEEPKS)
    list(APPEND objects
        LCAO_descriptor.cpp
        LCAO_descriptor_old.cpp
        LCAO_descriptor_dV.cpp
        LCAO_descriptor_io.cpp
        )
endif()

add_library(
    lcao
    OBJECT
    ${objects}
)

