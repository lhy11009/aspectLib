#!/bin/bash

################################################################################
# Tests functions from (base)
#
# Stdout:
#   number of successful tests and failed test to the terminal
#
# Run:
#   ./(test_base).sh
#
################################################################################

# source "${ASPECT_LAB_DIR}/bash_scripts/foo.sh"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"  # dir of this script
test_output_dir="${ASPECT_LAB_DIR}/.test"  # output to this dir
test_dir="${dir}/${FUNCNAME[0]}"


test_foo()
{
    ###
    # test function (foo)
    ###
    local_passed_tests=0
    local_failed_tests=0
   
    # compare_outputs "${FUNCNAME[0]}" "(standard outputs)" "(function outputs)"  # compare outputs
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message to terminal: numbers of successful and failed tests
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}


main(){
    ###
    # main function
    ###
    
    # parse
    # project=$1
    # server_info=$2

    # numbers of passed tests and failed tests
    passed_tests=0
    failed_tests=0
    
    # test parse_block_output_value
    test_parse_block_output
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # test parse_block_output_to_file
    # todo
    test_parse_block_output_to_file
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # message
    final_message 'test_utilities.sh' ${passed_tests} ${failed_tests}
}

main $@