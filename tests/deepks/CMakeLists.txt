enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME deepks_test_with_asan
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/deepks
    )

else()
    add_test(
        NAME deepks_test
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/deepks
    )
endif()
