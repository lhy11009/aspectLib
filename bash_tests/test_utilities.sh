#!/bin/bash

################################################################################
# Tests functions for aspect_lib.sh
# Run:
#   ./test_aspect_lib.sh
# Stdout:
#   test results
################################################################################

source "${ASPECT_LAB_DIR}/utilities.sh"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


################################################################################
# test bash_to_python_array()
# Inputs:
#   local_passed_tests: number of passed tests
#   local_failed_tests: number of failed tests
test_bash_to_python_array()
{
    local_passed_tests=0
    local_failed_tests=0
    
    # test 0
    bash_array=(0 100 200)
    bash_to_python_array
    compare_outputs "${FUNCNAME[0]}" "[0, 100, 200]" "${python_array_like}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}

################################################################################
# todo
# test python_to_bash_array
# Inputs:
#   local_passed_tests: number of passed tests
#   local_failed_tests: number of failed tests
test_python_to_bash_array()
{
    local_passed_tests=0
    local_failed_tests=0
    
    # test 0
    python_array_like="[0, 100, 200]"
    python_to_bash_array
    compare_outputs "${FUNCNAME[0]}" "0 100 200" "${bash_array[*]}"
    if [[ $? = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}


################################################################################
# todo
# test read_json_file
# read_json_file()
# Inputs:
#   local_passed_tests: number of passed tests
#   local_failed_tests: number of failed tests
test_read_json_file()
{
    local_passed_tests=0
    local_failed_tests=0

    # test 0
    local unmatched_output=0
    filein="./test1.json"
    # first key
    keys=("project")
    read_json_file
    compare_outputs "${FUNCNAME[0]}" "TwoDSubduction" "${value}"
    ((unmatched_output+=$?))
    # second key
    keys=("slurm" "nodes")
    read_json_file
    compare_outputs "${FUNCNAME[0]}" "1" "${value}"
    ((unmatched_output+=$?))
    # third key, this time an array
    keys=("slurm" "modules_to_unload")
    is_array='true'
    read_json_file
    compare_outputs "${FUNCNAME[0]}" "openmpi openmpi1" "${value[*]}"
    ((unmatched_output+=$?))
    if [[ ${unmatched_output} = 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}


################################################################################
# Test functions
test_element_in(){
	local _test_array=('a' 'b' 'c d')
	if ! element_in 'a' "${_test_array[@]}"; then
		cecho ${BAD} "test_element_in failed, 'a' is not in ${_test_array}[@]"
	fi
	if element_in 'c' "${_test_array[@]}"; then
		cecho ${BAD} "test_element_in failed, 'c' is in ${_test_array}[@]"
	fi
	cecho ${GOOD} "test_element_in passed"

}


test_parse_stdout(){
	# test the parse_stdout function, return values are last timestpe and time
	local _ifile="tests/integration/fixtures/task-2009375.stdout"
	if ! [[ -e ${_ifile} ]]; then
		cecho ${BAD} "test_parse_stdout failed, no input file ${_ifile}"
		exit 1
	fi
	parse_stdout ${_ifile}  # parse this file
	if ! [[ ${last_time_step} = "10" ]]; then
		cecho ${BAD} "test_parse_stdout failed, time_step is wrong"
		exit 1
	fi
	if ! [[ ${last_time} = "101705years" ]]; then
		cecho ${BAD} "test_parse_stdout failed, time is wrong"
		exit 1
	fi
	cecho ${GOOD} "test_parse_stdout passed"
}


test_read_log(){
	local log_file="${test_fixtures_dir}/test.log"
	read_log "${log_file}"
	if ! [[ "${return_value0}" = "tests/integration/fixtures tests/integration/fixtures" && "${return_value1}" = "2009375 2009376" ]]; then
		cecho ${BAD} "test_read_log failed, return values are not correct"
		return 1
	fi
	cecho ${GOOD} "test_read_log passed"
}


test_write_log(){
    local _ofile="${test_dir}/test.log"
    if [[ -e ${_ofile} ]]; then
        # remove older file
        eval "rm ${_ofile}"
    fi
    # test 1, write a non-existent job, it should return a NA status
    write_log_header "${_ofile}"
    write_log "${test_fixtures_dir}" "2009375" "${_ofile}"
    if ! [[ -e "${_ofile}" ]]; then
        cecho ${BAD} "test_write_log fails for test1, \"${_ofile}\"  doesn't exist"
	exit 1
    fi
    _output=$(cat "${_ofile}" | sed -n '2'p)
    if ! [[ ${_output} = "${test_fixtures_dir} 2009375 NA 10 101705years" ]]
    then
        cecho ${BAD} "test_write_log fails for test2, output format is wrong"
	exit 1
    fi
    cecho ${GOOD} "test_write_log passed"

}

test_clean_log(){
	local log_source_file="${test_fixtures_dir}/test.log"
	# copy source file to log file
	local log_file="${test_dir}/test.log"
	eval "cp ${log_source_file} ${log_file}"
	# call function
	clean_log "tests/integration/fixtures" "${log_file}"
	contents=$(cat "${log_file}")  # debug
	# compare file content
	if ! [[ "${contents}" = "job_dir job_id ST last_time_step last_time" ]]; then
		cecho ${BAD} "test_clean_log failed, file contents are not correct"
		return 1
	fi
	cecho ${GOOD} "test_clean_log passed"
}

test_fix_route() {
    local dir=$(pwd)
    # test1 test for relacing '~'
    fixed_route=$(fix_route "~/foo/ffoooo")
    [[ "${fixed_route}" = "${HOME}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 1"; exit 1; }
    # test2, test for replacing '.'
    fixed_route=$(fix_route "./foo/ffoooo")
    [[ "${fixed_route}" = "${dir}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 2"; exit 1; }
    # test3, test for replacing relative route
    fixed_route=$(fix_route "foo/ffoooo")
    [[ "${fixed_route}" = "${dir}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 3"; exit 1; }
    # test4, test for replacing '.'
    fixed_route=$(fix_route ".")
    [[ "${fixed_route}" = "${dir}" ]] || { cecho ${BAD} "test_fix_route failed for test 4"; exit 1; }
    # test5, test for relative address starts with .'
    fixed_route=$(fix_route ".test")
    [[ "${fixed_route}" = "${dir}/.test" ]] || { cecho ${BAD} "test_fix_route failed for test 5"; exit 1; }
    cecho ${GOOD} "test_fix_route passed"
}


################################################################################
# main function
# do all tests
# todo
main(){
    # parse
    project=$1
    server_info=$2

    passed_tests=0
    failed_tests=0

    # test python_to_bash_array
    test_python_to_bash_array
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # test bash_to_python_array
    test_bash_to_python_array
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # test read_json_file
    test_read_json_file
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # Test creating case and group with aspect_lib and submit to server
    if [[ -n $server_info ]]; then
        echo "0"
        # test_aspect_create_server
        # ((passed_tests+=local_passed_tests))
        # ((failed_tests+=local_failed_tests))
    else
        cecho ${WARN} "test_utilities.sh: no server info given, only do local tests"
    fi

    # message
    final_message 'test_utilities.sh' ${passed_tests} ${failed_tests}
}

main $@