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
# This group of tests test creating cases and groups under a project
# Cases are generated relative to the 'base.prm' file in the
# directory of "files/${project}". Cases are removed automatically
# after the tests succeed
test_aspect_create(){
    local_passed_tests=0
    local_failed_tests=0

    # test 0
    test_aspect_create0
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi

    # test1
    test_aspect_create1
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi

    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}
    return 0
}


################################################################################
# first test
# creat a case
test_aspect_create0(){
    local local_root=$(eval "echo \${${project}_DIR}")
    
    # case name and directory
    case_name="test0ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    
    # test file 
    source_dir="${dir}/test_aspect_lib/test_aspect_create"
    test_json_file="${source_dir}/config_case0.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    cp ${test_json_file} "${ASPECT_LAB_DIR}/config_case.json"
    
    # call function
    eval "${ASPECT_LAB_DIR}/aspect_lib.sh ${project} create"
    [[ $? -eq 0 ]] || { cecho ${BAD} "${FUNCNAME[0]}: create case failed with non-zero return value"; return 1; }

    # check case dir
    [[ -d ${case_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: create case failed. The directory of created case doesn't exist"; return 1; }

    # check .prm file
    prm_file="${case_dir}/case.prm"
    prm_file_std="${source_dir}/case0.prm"
    compare_files "${FUNCNAME[0]}" "${prm_file_std}" "${prm_file}"
    [[ $? -eq 0 ]] || return 1
    
    return 0
}


################################################################################
# second test
# creat a group
test_aspect_create1(){
    local local_root=$(eval "echo \${${project}_DIR}")
    
    # group name and directory
    group_name='test_group'
    group_dir="${local_root}/${group_name}"
    [[ -d "${group_dir}" ]] && eval "rm -r ${group_dir}"  # remove older dir

    # test_file
    source_dir="${dir}/test_aspect_lib/test_aspect_create"
    test_json_file="${source_dir}/config_group1.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    cp ${test_json_file} "${ASPECT_LAB_DIR}/config_group.json"

    # call function
    eval "${ASPECT_LAB_DIR}/aspect_lib.sh ${project} create_group"
    [[ $? -eq 0 ]] || { cecho ${BAD} "${FUNCNAME[0]}: create group failed with non-zero return value"; return 1; }

    # check group dir
    [[ -d ${group_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: create group failed. The directory of created group doesn't exist"; return 1; }

    # check sub dirs
    sub_dirs=("test1ULV3.000e+01testIAR6" "test1ULV1.000e+02testIAR6" "test1ULV3.000e+01testIAR8" "test1ULV1.000e+02testIAR8")
    for sub_dir in ${sub_dirs[@]}; do
        [[ -d "${group_dir}/${sub_dir}" ]] || { cecho ${BAD} "${FUNCNAME[0]}: create group failed. The subdirectory of created group: ${sub_dir} doesn't exist"; return 1; }
    done
    
    # check one of .prm files
    prm_file="${group_dir}/test1ULV3.000e+01testIAR6/case.prm"
    prm_file_std="${source_dir}/case1.prm"
    compare_files "${FUNCNAME[0]}" "${prm_file_std}" "${prm_file}"
    [[ $? -eq 0 ]] || return 1
    
    return 0
}


################################################################################
# This group of tests test creating cases and groups under a project and then submit to server
# Cases are generated relative to the 'base.prm' file in the
# directory of "files/${project}". Cases are removed automatically
# after the tests succeed
test_aspect_create_server(){
    local_passed_tests=0
    local_failed_tests=0

#    # test 0
#    test_aspect_create_server0
#    if [[ $? -eq 0 ]]; then
#        ((local_passed_tests++))
#    else
#        ((local_failed_tests++))
#    fi
    
    # test 1
    test_aspect_create_server1
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi

    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}

    return 0
}


################################################################################
# first test
# creat a case and submit to server
test_aspect_create_server0(){
    local local_root=$(eval "echo \${${project}_DIR}")
    
    # case name and directory
    case_name="test0ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    
    # test file 
    source_dir="${dir}/test_aspect_lib/test_aspect_create"
    test_json_file="${source_dir}/config_case0.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    cp ${test_json_file} "${ASPECT_LAB_DIR}/config_case.json"

    # call function
    eval "${ASPECT_LAB_DIR}/aspect_lib.sh ${project} create_submit ${server_info}"
    [[ $? -eq 0 ]] || { cecho ${BAD} "${FUNCNAME[0]}: create case failed with non-zero return value"; return 1; }

    # get job id
    local job_id=$(cat ".temp")
    [[ ${job_id} =~ ^[0-9]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: submit case failed with invalid job id"; return 1; }

    # cancel job and delete case folder when done
    ssh "${server_info}" eval "\${SCANCEL} ${job_id}"
    
    return 0
}


################################################################################
# second test
# creat a case and submit to server with a log file generated
test_aspect_create_server1(){
    local local_root=$(eval "echo \${${project}_DIR}")
    local_log_file="${ASPECT_LAB_DIR}/.test/test.log"
    
    # case name and directory
    case_name="test0ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    
    # test file 
    source_dir="${dir}/test_aspect_lib/test_aspect_create"
    test_json_file="${source_dir}/config_case0.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    cp ${test_json_file} "${ASPECT_LAB_DIR}/config_case.json"
    
    # todo clean previous files
    eval "${ASPECT_LAB_DIR}/process.sh remove ${local_log_file} ${server_info}"

    # call function
    eval "${ASPECT_LAB_DIR}/aspect_lib.sh ${project} create_submit ${server_info} ${local_log_file}"
    [[ $? -eq 0 ]] || { cecho ${BAD} "${FUNCNAME[0]}: create case failed with non-zero return value"; return 1; }

    # get job id
    local job_id=$(cat ".temp")
    [[ ${job_id} =~ ^[0-9]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: submit case failed with invalid job id"; return 1; }

    # cancel job and delete case folder when done
    ssh "${server_info}" eval "\${SCANCEL} ${job_id}"
    
    # scp log file
    eval "${ASPECT_LAB_DIR}/process.sh update_from_server ${local_log_file} ${server_info}"
    
    # check log file
    [[ -e "${local_log_file}" ]] || { cecho ${BAD} "${FUNCNAME[0]}: submit case failed: log file is not generated"; return 1; }

    # check content of log file 
    # get second line
    local line2=$(sed -n '2'p "${local_log_file}")
    [[ "${line2}" =~ ^.*test0ULV1\.000e\+02testIAR8\ [0-9]* ]] || { cecho ${BAD} "${FUNCNAME[0]}: submit case failed: content of log file is not correct"; return 1; }
    
    return 0
}


# This test test creating cases and groups under a project
# Cases are generated relative to the 'base.prm' file in the
# directory of "files/${project}".
# Cases are then submited to a server to run.
# Cases are removed automatically, after the tests succeed
test_aspect_lib_remote(){


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
# Test translating a visit script
test_translate_visit()
{
    local_passed_tests=0
    local_failed_tests=0

    # test 0
    test_translate_visit0
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi

    # message
    final_message ${FUNCNAME[0]} ${local_passed_tests} ${local_failed_tests}

}


################################################################################
# Test translating a visit script: test 0
test_translate_visit0()
{
    # todo
    # A standard interface to do integration
    source_dir="${dir}/test_aspect_lib/test_translate_visit"

    # prepare input file
    filein="${source_dir}/visit_keys_values"
    filein1="${source_dir}/initial_slab.py"
    cp "${filein}"  "$ASPECT_LAB_DIR/"

    # output file path
    file_out_std="${source_dir}/initial_slab_std.py"
    file_out="${ASPECT_LAB_DIR}/visit_scripts_temp/initial_slab.py"

    # run command
    eval "${ASPECT_LAB_DIR}/aspect_lib.sh ${project} translate_visit $filein1"

    compare_files "${FUNCNAME[0]}" "${file_out_std}" "${file_out}"
    [[ $? -eq 0 ]] || return 1
    return 0
}


################################################################################
# main function
# do all tests
main(){
    # parse
    project=$1
    server_info=$2

    passed_tests=0
    failed_tests=0

    # Test translating visit scripts
    test_translate_visit
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # Test creating case and group with aspect_lib
    test_aspect_create
    ((passed_tests+=local_passed_tests))
    ((failed_tests+=local_failed_tests))

    # Test creating case and group with aspect_lib and submit to server
    if [[ -n $server_info ]]; then
        test_aspect_create_server
        ((passed_tests+=local_passed_tests))
        ((failed_tests+=local_failed_tests))
    else
        cecho ${WARN} "test_aspect_lib.sh: no server info given, only do local tests"
    fi

    # message
    final_message 'test_aspect_lib.sh' ${passed_tests} ${failed_tests}
}

main $@