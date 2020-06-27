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
    # submit using slurm.sh,
    # determine if there is a valid job id, todo
    ssh ${server_info} << EOF > '.temp'
        eval "slurm.sh -n ${total_tasks} -t ${time_by_hour} -p ${partition} ${remote_case_prm}"
EOF
    # get job_id
    local _info=$(cat '.temp'| sed -n '$'p)
    local job_id=$(echo "${_info}" | sed 's/Submitted\ batch\ job\ //')
    if ! [[ ${job_id} != '' && ${job_id} =~ ^[0-9]*$  ]]; then
        cecho ${BAD} "submit case: $(case_name) failed"
        return 1
    else
        cecho ${GOOD} "submit case: ${case_name} succeeded, job id: ${job_id}"
        return 0
    fi
}


test(){
    # do test, todo
    # create jobs
    eval "./aspect_lib.sh TwoDSubduction create"
    # create a group of jobs
    eval "./aspect_lib.sh TwoDSubduction create_group"
    # submit jobs
    eval "./aspect_lib.sh TwoDSubduction submit ULV1.000e+02testIAR8 lochy@peloton.cse.ucdavis.edu"
}


main(){
    # parameter list, todo
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local py_script="shilofue.${project}"
    # check project
    [[ -d ${local_root} ]] || { cecho ${BAD} "Project ${project} is not supported"; exit 1; }
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
        submit "${case_dir}" "${remote_case_dir}" "${server_info}"
    elif [[ ${_commend} = 'submit_group' ]]; then
        # todo
        local group_name="$3"
        local server_info="$4"
        local group_dir="${local_root}/${case_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        # get a list of cases, todo
        # local case_dirs=
        for case_dir in ${case_dirs[@]}; do
            local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"}
            # call submit functions
            submit "${case_dir}" "${remote_case_dir}" "${server_info}"
        done
    elif [[ ${_commend} = 'create_submit' ]]; then
        # todo
        # get remote foler
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"}
        # submit to server, todo
        ssh ${server_info} << EOF
            eval "slurm"
EOF
    elif [[ ${_commend} = 'create_submit_group' ]]; then
        # todo
        echo "foo"
    elif [[ ${_commend} = 'test' ]]; then
        # do test, todo
        echo "foo"
    else
        cecho ${BAD} "Bad commend: ${_commend}"
    fi
    return 0
}

main $@
