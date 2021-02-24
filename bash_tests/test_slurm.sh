#!/bin/bash

################################################################################
# Tests functions for aspect_lib.sh
# Run:
#   ./test_slurm.sh TwoDSubduction lochy@peloton.cse.ucdavis.edu
# Stdout:
#   test results
################################################################################

source "${ASPECT_LAB_DIR}/utilities.sh"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

################################################################################
# first test
# creat a case and submit to server with a log file generated
test_slurm_submit_server1(){
    source_dir="${dir}/test_slurm/slurm_test_submit1"
    [[ -d ${source_dir} ]] || cecho $BAD "source_dir doesn't exist"
    # todo
    local project_dir="${ASPECT_PROJECT_DIR}/${project}"
    # get remote variables
    get_remote_environment "${server_info}" "${project}_DIR"
    local remote_root=${return_value}
    get_remote_environment "${server_info}" "ASPECT_LAB_DIR"
    local remote_lib_dir=${return_value}
    
    # make case P16
    local remote_source_dir="${remote_lib_dir}/bash_tests/test_slurm/slurm_test_submit1"
    local remote_case_dir="${remote_root}/slurm_test_submit1" # substitution
    local remote_case_prm="${remote_case_dir}/case.prm"
    local remote_log_file="${remote_lib_dir}/.output/job.sh"
    local remote_time_file="${remote_lib_dir}/machine_time"
    
    # generate job.sh file
    local addition=""
    local flag=""
    ssh ${server_info} << EOF > ".temp"
        eval "[[ -d ${remote_case_dir} ]] && rm -r ${remote_case_dir}"
        eval "cp -r ${remote_source_dir} ${remote_case_dir}"
        eval "slurm.sh -N 1 -n 16 -t 24 ${addition} -l ${remote_log_file} -lt ${remote_time_file} -P ${project} ${flag} ${remote_case_prm}"
EOF
    local _content=$(cat ".temp")
    [[ ${_content} =~ "Failure" ]] && { cecho "${BAD}" "${FUNCNAME[0]}, failure in ssh process"; return 1; }
    
    # get job id
    local _info=$(cat '.temp'| sed -n '$'p)
    local job_id=$(echo "${_info}" | sed 's/Submitted\ batch\ job\ //')
    [[ ${job_id} =~ ^[0-9]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: submit case failed with invalid job id"; return 1; }

    # check time file  
    ssh ${server_info} << EOF > ".temp"
        eval " [[ -e ${remote_time_file} ]] && echo \"${remote_time_file} generated successfully\" "
EOF
    local _info=$(cat '.temp'| sed -n '$'p)
    [[ ${_info} =~ "successfully" ]] || { cecho ${BAD} "${FUNCNAME[0]}: remote time file is not generated"; return 1; }

    # cancel job and delete case folder when done
    ssh "${server_info}" eval "\${SCANCEL} ${job_id}"
}


################################################################################
# second test
# creat a case and submit to server with a log file generated
# with a hold option
test_slurm_submit_server2(){
    source_dir="${dir}/test_slurm/slurm_test_submit2"
    [[ -d ${source_dir} ]] || cecho $BAD "source_dir doesn't exist"
    # todo
    local project_dir="${ASPECT_PROJECT_DIR}/${project}"
    # get remote variables
    get_remote_environment "${server_info}" "${project}_DIR"
    local remote_root=${return_value}
    get_remote_environment "${server_info}" "ASPECT_LAB_DIR"
    local remote_lib_dir=${return_value}
    
    # make case P16
    local remote_source_dir="${remote_lib_dir}/bash_tests/test_slurm/slurm_test_submit2"
    local remote_case_dir="${remote_root}/slurm_test_submit2" # substitution
    local remote_case_prm="${remote_case_dir}/case.prm"
    local remote_log_file="${remote_lib_dir}/.output/job.sh"
    local remote_time_file="${remote_lib_dir}/machine_time"
    
    # generate job.sh file
    local addition=""
    local flag="--hold"  # add a hold option
    ssh ${server_info} << EOF > ".temp"
        eval "[[ -d ${remote_case_dir} ]] && rm -r ${remote_case_dir}"
        eval "cp -r ${remote_source_dir} ${remote_case_dir}"
        eval "slurm.sh -N 1 -n 16 -t 24 ${addition} -l ${remote_log_file} -lt ${remote_time_file} -P ${project} ${flag} ${remote_case_prm}"
EOF
    local _content=$(cat ".temp")
    [[ ${_content} =~ "Failure" ]] && { cecho "${BAD}" "${FUNCNAME[0]}, failure in ssh process"; return 1; }

    # scp the job.sh file
    local job_sh_file="${ASPECT_LAB_DIR}/.test/job.sh"
    eval "${RSYNC} ${server_info}:${remote_case_dir}/job.sh ${job_sh_file}"

    [[ "${server_info}" =~ "peloton" ]] && job_sh_file_std="${source_dir}/peloton_job_std.sh"
    [[ -e ${job_sh_file_std} ]] || { cecho "${BAD}" "${FUNCNAME[0]}, job_sh_file_std doesn't exist."; return 1; }
    compare_files "${FUNCNAME[0]}" ${job_sh_file_std} ${job_sh_file}
}


################################################################################
# This group of tests submit to server with slurm.sh
# Cases are removed automatically
# after the tests succeed
test_slurm_submit_server(){
    local_passed_tests=0
    local_failed_tests=0
    
    # test 1
    test_slurm_submit_server1
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi
    
    # test 2
    test_slurm_submit_server2
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
# main function
# do all tests
main(){
    # parse
    project=$1
    set_server_info "$2"

    passed_tests=0
    failed_tests=0

    # sub-environment variable 
    local local_root=$(eval "echo \${${project}_DIR}")

    # Test creating case and group with aspect_lib and submit to server
    if [[ -n $server_info ]]; then
        test_slurm_submit_server
        ((passed_tests+=local_passed_tests))
        ((failed_tests+=local_failed_tests))
    else
        cecho ${BAD} "test_slurm.sh: no server info given"
    fi

    # message
    final_message 'test_slurm.sh' ${passed_tests} ${failed_tests}
}

main $@