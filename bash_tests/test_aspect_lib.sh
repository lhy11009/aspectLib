#!/bin/bash

################################################################################
# Tests functions
################################################################################


# This test test creating cases and groups under a project
# Cases are generated relative to the 'base.prm' file in the
# directory of "files/${project}". Cases are removed automatically
# after the tests succeed
test_aspect_lib(){
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")

    # test1 create a case
    case_name="baseULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create'
    quit_if_fail "test_aspect_lib: create case failed"
    eval "rm -r ${case_dir}" # remove dir after test suceed
    
    # test2 create a group
    group_name='test_group'
    group_dir="${local_root}/${group_name}"
    [[ -d "${group_dir}" ]] && eval "rm -r ${group_dir}"  # remove older dir
    # test_file
    test_json_file="${dir}/tests/integration/fixtures/config_group.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create_group'
    quit_if_fail "test_aspect_lib: create group failed"
    cecho ${GOOD} "test_aspect_lib succeeded"
    eval "rm -r ${group_dir}" # remove dir after test succeeded
}


# This test test creating cases and groups under a project
# Cases are generated relative to the 'base.prm' file in the
# directory of "files/${project}".
# Cases are then submited to a server to run.
# Cases are removed automatically, after the tests succeed
test_aspect_lib_remote(){
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local server_info="$2"

    # test1 create and submit case
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    case_name="baseULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    # ./process.sh clean "${case_dir}" "${server_info}"
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create_submit' "${server_info}"
    quit_if_fail "test_aspect_lib_remote submit case failed"
    local job_id=$(cat ".temp")
    # cancel job and delete case folder when done
    ssh "${server_info}" eval "\${SCANCEL} ${job_id}"
    eval "rm -r ${case_dir}"

    # test2 create and submit case with a log file
    case_name="baseULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    local_log_file=".test/test.log"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # todo clean previous files
    ./process.sh remove "${local_log_file}" "${server_info}"
    # call function
    ./aspect_lib.sh "${project}" 'create_submit' "${server_info}" "${local_log_file}"
    quit_if_fail "test_aspect_lib_remote submit case with log file failed"
    local job_id=$(cat ".temp")
    # cancel job and delete case folder when done
    ssh "${server_info}" eval "\${SCANCEL} ${job_id}"
    eval "rm -r ${case_dir}"
    # scp log file
    ./process.sh update_from_server "${local_log_file}" "${server_info}"
    # check log file
    if ! [[ -e "${local_log_file}" ]]; then
        cecho ${BAD} "test_aspect_lib_remote submit case with log file failed log_file ${local_log_file} is not generated"
        exit 1
    fi
    local line2=$(sed -n '2'p "${local_log_file}")  # get second line of the file
    if ! [[ "${line2}" =~ ^.*baseULV1\.000e\+02testIAR8\ [0-9]* ]]; then
        cecho ${BAD} "test_aspect_lib_remote submit case with log file failed, content in log_file ${local_log_file} is not correct"
        exit 1
    fi

    # test3 create and submit group
    group_name='test_group'
    group_dir="${local_root}/${group_name}"
    [[ -d "${group_dir}" ]] && eval "rm -r ${group_dir}"  # remove older dir
    # test_file
    test_json_file="${dir}/tests/integration/fixtures/config_group.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create_submit_group' "${server_info}"
    quit_if_fail "test_aspect_lib_remote submit group failed"
    # get job_id and cancel cases
    local job_ids
    local job_id
    IFS=' ' read -r -a job_ids <<< $(cat ".temp")
    for job_id in ${job_ids[@]}; do
        ssh "${server_info}" eval "\${SCANCEL} ${job_id}"  # cancel job after test
    done
    cecho ${GOOD} "test_aspect_lib_remote succeeded"
}

################################################################################
# Compare file content, return 1 if the contents are different
# inputs:
#    $1 : file1
#    $2 : file2
# stderr:
#   if two files are different
# return:
#   0 : succeed
#   1 : failed
# output:
#   message:
#       message to give when test failed, other wise ''
check_file_content(){
    message=""

    # todo
    file1="$1"
    file2="$2"
    file_difference = $(diff "${file1}" "${file2}")
    if ! [[ -z ${file_difference} ]]; then
        message="${message} ${file1} and ${file2} are different"
        cecho ${BAD} "${message} Difference are:"
        cecho ${BAD} "${message} ${file_difference}"
        return 1
    fi
}


################################################################################
# Test translating a visit script
test_translate_visit()
{
    # todo
    # A standard interface to do integration
    local test_passed=0
    local test_failed=0
    local error_message=""

    # test 1
    check_file_content "${file_out_std}" "${file_out}"
    if [ $? -eq 0 ]; then
        ((++test_passed))
    else
        ((++test_failed))
        error_message="${error_messge}, ${message}"
    fi


    cecho ${GOOD} "test_translate_visit: ${tests_passed} tests passed"
    
    if ! [[ ${test_failed} -eq 0 ]]; then
        cecho ${BAD} "test_translate_visit: ${test_failed} tests failed"
        cecho ${BAD} "Error message is"
        cecho ${BAD} "${error_message}" 
    fi
    
    echo '0'
}

main(){
    test_translate_visit()
}
