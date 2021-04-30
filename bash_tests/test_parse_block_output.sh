#!/bin/bash

################################################################################
# Tests functions for parse_block_output.sh
# Run:
#   ./test_parse_block_output.sh
# Stdout:
#   test results
################################################################################

source "${ASPECT_LAB_DIR}/utilities.sh"
source "${ASPECT_LAB_DIR}/bash_scripts/parse_block_output.sh"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
test_output_dir="${ASPECT_LAB_DIR}/.test"
test_dir="${dir}/test_parse_block_output"


#################################################################################
# test parse_block_output()
# Inputs:
test_parse_block_output()
{
    local_passed_tests=0
    local_failed_tests=0
   
    # test 1: parse the block output at the end of a log file
    local log_file="${test_dir}/log.txt"
    parse_block_output_wallclock "${log_file}"
    compare_outputs "${FUNCNAME[0]}" "497s 950s 950s" "${return_values[*]}"
    parse_block_output "${log_file}" "Assemble Stokes system" 
    compare_outputs "${FUNCNAME[0]}" "155s 314s 314s" "${return_values[*]}"
    parse_block_output "${log_file}" "Initialization" 
    compare_outputs "${FUNCNAME[0]}" "0.892s 0.892s 0.892s" "${return_values[*]}"
    parse_block_output "${log_file}" "Setup matrices" 
    compare_outputs "${FUNCNAME[0]}" "18.2s 18.2s 18.2s" "${return_values[*]}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # test 2: parse the block output at the end of a log file, contains a e+00 input this time
    local log_file="${test_dir}/log1.txt"
    parse_block_output_wallclock "${log_file}"
    compare_outputs "${FUNCNAME[0]}" "534s 1.05e+03s 1.05e+03s" "${return_values[*]}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}


#################################################################################
# test parse_block_output_to_file()
# todo
# Inputs:
test_parse_block_output_to_file()
{
    local_passed_tests=0
    local_failed_tests=0

    # test 1, parse a file with two keys 
    local log_file="${test_dir}/log.txt"
    local ofile="${test_output_dir}/block_output.txt"
    parse_block_output_to_file "${log_file}" "${ofile}" "Assemble Stokes system" "Setup matrices"
    # compare
    std_ofile="${test_dir}/block_output_std.txt"
    [[ -e ${std_ofile} ]] || cecho ${BAD} "${FUNC_NAME[0]}: std_ofile doesn't exists"
    compare_files "${FUNCNAME[0]}" "${std_ofile}" "${ofile}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # test 2, parse a file with two keys, this time with entry like e+01
    local log_file="${test_dir}/log1.txt"
    local ofile="${test_output_dir}/block_output1.txt"
    parse_block_output_to_file "${log_file}" "${ofile}" "Assemble Stokes system" "Setup matrices"
    # compare
    std_ofile="${test_dir}/block_output_std1.txt"
    [[ -e ${std_ofile} ]] || cecho ${BAD} "${FUNC_NAME[0]}: std_ofile doesn't exists"
    compare_files "${FUNCNAME[0]}" "${std_ofile}" "${ofile}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # test 3, parse a file, looking for steps and times
    local log_file="${test_dir}/log2.txt"
    local ofile="${test_output_dir}/block_output2.txt"
    eval "awk -f ${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output ${log_file} > ${ofile}"
    # compare
    std_ofile="${test_dir}/block_output_std2.txt"
    [[ -e ${std_ofile} ]] || cecho ${BAD} "${FUNC_NAME[0]}: std_ofile doesn't exists"
    compare_files "${FUNCNAME[0]}" "${std_ofile}" "${ofile}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
    
}

main(){
    # parse
    project=$1
    server_info=$2

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
