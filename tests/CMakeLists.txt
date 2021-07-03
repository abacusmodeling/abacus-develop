enable_testing()

find_program(BASH bash)

add_test(
    NAME integrated_test
    COMMAND ${BASH} Autotest.sh
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}
)
