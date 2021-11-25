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
test_dir="${dir}/test_parse_case"

test_parse_block_output()
{
    ###
    # test parse_block_output with the awk file parse_block_output
    # this make use of the block output of aspect (instead of the time steps)
    # Results shown are "step, time, wall clock"
    ###
    local_passed_tests=0
    local_failed_tests=0
    # test 1: this is file that contains block output, check the values of outputs
    local log_file="${test_dir}/log.txt"
    [[ -e "${log_file}" ]] || { cecho "$BAD" "No such file ${log_file}"; exit 1; }
    local temp=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output" "${log_file}")
    [[ $? = 0 ]] ||  cecho "$BAD" "The previous command doesn't work."
    IFS=$'\n'; local entries=(${temp})
    local outputs="${entries[*]: -1}"
    [[ ${outputs} =~ "#" ]] && outputs=""  # if only header is retured, reset
    compare_outputs "${FUNCNAME[0]}" "10      70508.1   2.68e+03" "${outputs}"  # compare outputs
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    # test 2: this is file that doesn't contain block output, check that output is vacant
    local log_file="${test_dir}/log1.txt"
    [[ -e "${log_file}" ]] || { cecho "$BAD" "No such file ${log_file}"; exit 1; }
    local temp=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output" "${log_file}")
    [[ $? = 0 ]] ||  cecho "$BAD" "The previous command doesn't work."
    IFS=$'\n'; local entries=(${temp})
    local outputs="${entries[*]: -1}"
    [[ ${outputs} =~ "#" ]] && outputs=""  # if only header is retured, reset
    compare_outputs "${FUNCNAME[0]}" "" "${outputs}"  # compare outputs
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    # message to terminal: numbers of successful and failed tests
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}

test_parse_timestep_last()
{
    ###
    # test with the awk file parse_log_last_step
    # this make use of the timestep output of aspect.
    # Results shown are "step, time"
    ###
    local_passed_tests=0
    local_failed_tests=0
    # test 1
    local log_file="${test_dir}/log.txt"
    [[ -e "${log_file}" ]] || { cecho "$BAD" "No such file ${log_file}"; exit 1; }
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_log_last_step" "${log_file}")
    [[ $? = 0 ]] ||  cecho "$BAD" "The previous command doesn't work."
    compare_outputs "${FUNCNAME[0]}" "10      70508.1" "${outputs}"  # compare outputs
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    # test 2, this one doesn't contain block output (once every 10 steps) but contains time step
    # outputs nonetheless. This intends to show the different usage of these two.
    local log_file="${test_dir}/log1.txt"
    [[ -e "${log_file}" ]] || { cecho "$BAD" "No such file ${log_file}"; exit 1; }
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_log_last_step" "${log_file}")
    [[ $? = 0 ]] ||  cecho "$BAD" "The previous command doesn't work."
    compare_outputs "${FUNCNAME[0]}" "1      9791.84" "${outputs}"  # compare outputs
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}

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

    # test parse_timestep_last
    test_parse_timestep_last
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))


    # message
    final_message 'test_parse_case.sh' ${passed_tests} ${failed_tests}
}

main $@