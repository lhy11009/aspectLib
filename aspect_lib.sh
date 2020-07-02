#!/bin/bash
# case manager
# Usage:
#   ./aspect_lib.sh + command + options

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"

source "${dir}/utilities.sh"
source "${dir}/extra_configurations.sh"

create_case(){
    # create a case locally
    local py_script="$1"
    local local_root="$2"
    eval "python -m ${py_script} create -j config_case.json 2>&1 > .temp"
    [[ $? -eq 1 ]] && { cecho ${BAD} "${py_script} failed"; exit 1; }
    # get case name
    _info=$(cat ".temp")
    case_name=$(echo "${_info}" | sed -n '2'p)
    # assertion
    local case_dir="${local_root}/${case_name}"
    local case_prm="${case_dir}/case.prm"
    [[ -d ${case_dir} && -e ${case_prm} ]] || { cecho ${BAD} "Case generation failed"; exit 1; }
    cecho ${GOOD} "${_info}"
}


create_group(){
    # create a group of case locally
    local py_script="$1"
    local local_root="$2"
    eval "python -m ${py_script} create_group -j config_group.json 2>&1 > .temp"
    [[ $? -eq 1 ]] && { cecho ${BAD} "${py_script} failed"; exit 1; }
    # get case names
    _info=$(cat ".temp")
    local group_name=$(echo "${_info}" | sed -n '2'p)
    IFS=' ' read -r -a case_names <<< $(cat ".temp" | sed -n '1,2!'p)
    group_dir="${local_root}/${group_name}"
    for case_name in ${case_names[@]}; do
        local case_dir="${group_dir}/${case_name}"
        local case_prm="${case_dir}/case.prm"
        [[ -d ${case_dir} && -e ${case_prm} ]] || { cecho ${BAD} "Case generation failed"; exit 1; }
    done
    cecho ${GOOD} "${_info}"
}


submit(){
    local case_dir="$1"
    local case_name=$(basename "${case_dir}")
    local romote_case_dir="$2"
    local server_info="$3"
    local flag=''  # a vacant flag for adding optional parameters
    [[ ${RSYNC} = ''  ]] && local RSYNC='rsync'  # replace rsync with rsync1
    local case_prm="${case_dir}/case.prm"
    local remote_case_prm="${remote_case_dir}/case.prm"
    # get configuration from a file
    total_tasks=$(sed -n '1'p "slurm_config")
    time_by_hour=$(sed -n '2'p "slurm_config")
    partition=$(sed -n '3'p "slurm_config")
    # rsync to remote
    local remote_target=$(dirname "${remote_case_dir}")
    eval "${RSYNC} -avur ${case_dir} ${server_info}:${remote_target}"
    # add an optional log file
    [[ "$4" != '' ]] && flag="${flag} -l $4"  # add -l log_file to flag, if $4 given
    # submit using slurm.sh,
    # determine if there is a valid job id, todo
    ssh ${server_info} << EOF > '.temp'
        eval "slurm.sh -n ${total_tasks} -t ${time_by_hour} -p ${partition} ${remote_case_prm} ${flag}"
EOF
    # get job_id
    local _info=$(cat '.temp'| sed -n '$'p)
    local job_id=$(echo "${_info}" | sed 's/Submitted\ batch\ job\ //')
    if ! [[ ${job_id} != '' && ${job_id} =~ ^[0-9]*$  ]]; then
        cecho ${BAD} "submit case: ${case_name} failed"
        return 1
    else
        cecho ${GOOD} "submit case: ${case_name} succeeded, job id: ${job_id}"
        return 0
    fi
}


test_aspect_lib(){
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    # test1 create a case #############################################
    case_name="ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create'
    quit_if_fail "test_aspect_lib: create case failed"
    # test2 create a group ############################################
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
}

test_aspect_lib_remote(){
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local server_info="$2"
    # test1 create and submit case ####################################
    case_name="ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # call function
    ./aspect_lib.sh "${project}" 'create_submit' "${server_info}"
    quit_if_fail "test_aspect_lib_remote submit case failed"
    # test2 create and submit case with a log file######################
    case_name="ULV1.000e+02testIAR8"
    case_dir="${local_root}/${case_name}"
    local_log_file=".test/test.log"
    [[ -d "${case_dir}" ]] && eval "rm -r ${case_dir}"  # remove older dir
    # test file
    test_json_file="${dir}/tests/integration/fixtures/config_case.json"
    [[ -e "${test_json_file}" ]] || { cecho ${BAD} "test file ${test_json_file} doesn't exist"; exit 1; }
    eval "cp ${test_json_file} ${dir}/"
    # todo clean previous files
    ./process.sh clean "${local_log_file}" "${server_info}"
    # call function
    ./aspect_lib.sh "${project}" 'create_submit' "${server_info}" "${local_log_file}"
    quit_if_fail "test_aspect_lib_remote submit case with log file failed"
    # scp log file
    ./process.sh update_from_server "${local_log_file}" "${server_info}"
    # check log file
    if ! [[ -e "${local_log_file}" ]]; then
        cecho ${BAD} "test_aspect_lib_remote submit case with log file failed log_file ${local_log_file} is not generated"
        exit 1
    fi
    local line2=$(sed -n '2'p "${local_log_file}")  # get second line of the file
    if ! [[ "${line2}" =~ ^.*ULV1\.000e\+02testIAR8\ [0-9]* ]]; then
        cecho ${BAD} "test_aspect_lib_remote submit case with log file failed, content in log_file ${local_log_file} is not correct"
        exit 1
    fi
    # test3 create and submit group ####################################
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
    cecho ${GOOD} "test_aspect_lib_remote succeeded"
}


main(){
    # parameter list, todo
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local py_script="shilofue.${project}"
    # check project
    [[ -d ${local_root} ]] || { cecho ${BAD} "Project ${project} is not included"; exit 1; }
    # get commend
    _commend="$2"
    if [[ ${_commend} = 'create' ]]; then
        create_case "${py_script}" "${local_root}"
    elif [[ ${_commend} = 'create_group' ]]; then
        create_group "${py_script}" "${local_root}"
    elif [[ ${_commend} = 'submit' ]]; then
        local case_name="$3"
        local server_info="$4"
        local case_dir="${local_root}/${case_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"} # substitution
        local log_file="$5"  # add an optional log_file
        if [[ "${log_file}" != '' ]]; then
            # if there is no $5 given, log file is ''
            log_file=$(fix_route "${log_file}")
            log_file=${log_file/"${local_root}"/"${remote_root}"} # substitution
        fi
        submit "${case_dir}" "${remote_case_dir}" "${server_info}" "${log_file}"
    elif [[ ${_commend} = 'submit_group' ]]; then
        # todo
        local group_name="$3"
        local server_info="$4"
        local group_dir="${local_root}/${group_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        # get a list of cases and submit
        for case_dir in "${group_dir}/"*/; do
            local _files=$(ls "${case_dir}")
            if [[ "${_files[@]}" =~ 'case.prm' ]]; then
                local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"}
                # call submit functions
                submit "${case_dir}" "${remote_case_dir}" "${server_info}"
            fi
        done
        return 0
    elif [[ ${_commend} = 'create_submit' ]]; then
        local server_info="$3"
        local log_file="$4"  # optional log file
        ./aspect_lib.sh "${project}" 'create'
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh create failed"; exit 1; }
        # get case name
        local _info=$(cat ".temp")
        local case_name=$(echo "${_info}" | sed -n '2'p)
        # submit to server
        ./aspect_lib.sh "${project}" 'submit' "${case_name}" "${server_info}" "${log_file}"
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh submit failed"; exit 1; }
    elif [[ ${_commend} = 'create_submit_group' ]]; then
        local server_info="$3"
        ./aspect_lib.sh "${project}" 'create_group'
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh create failed"; exit 1; }
        # get group name
        local _info=$(cat ".temp")
        local group_name=$(echo "${_info}" | sed -n '2'p)
        # call self
        ./aspect_lib.sh "${project}" 'submit_group' "${group_name}" "${server_info}"
    elif [[ ${_commend} = 'terminate' ]]; then
        # todo
        echo '0'
    elif [[ ${_commend} = 'remove' ]]; then
        # todo
        echo '0'
    elif [[ ${_commend} = 'test' ]]; then
        # do test
        test_aspect_lib "${project}"
    elif [[ ${_commend} = 'remote_test' ]]; then
        # do test
        [[ "$#" -eq 3 ]] || { cecho ${BAD} "for remote_test, server_info must be given"; exit 1; }
        local server_info="$3"
        test_aspect_lib_remote "${project}" "${server_info}"
    else
        cecho ${BAD} "Bad commend: ${_commend}"
    fi
    return 0
}

main $@
