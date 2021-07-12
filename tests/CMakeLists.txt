enable_testing()

find_program(BASH bash)

add_test(
    NAME integrated_test
    COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 16
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}
)
