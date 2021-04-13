#!/bin/bash


################################################################################
# do affinity test
#
# Dependencies:
#    env:
#        ASPECT_LAB_DIR
#        ASPECT_SOURCE_DIR
#
# Example Usage:
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
# source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"

do_affinity_test(){
    ###
    # affinty test on server
    # 
    # name of tests are P*B#, * could be 16, 32, 64, 128
    # could be 1 or 2, this is wheter the binded option is used.
    ###
    source_dir="${ASPECT_LAB_DIR}/files/${project}/affinity_test"
    [[ -d ${source_dir} ]] || cecho $BAD "source_dir doesn't exist"
    # todo
    local project_dir="${ASPECT_PROJECT_DIR}/${project}"
    # get remote variables
    get_remote_environment "${server_info}" "${project}_DIR"
    local remote_root=${return_value}
    get_remote_environment "${server_info}" "ASPECT_LAB_DIR"
    local remote_lib_dir=${return_value}
    
    if [[ "${server_info}" =~ "peloton" ]]; then
        target_dir="${project_dir}/peloton_affinity_test"
    else
        target_dir="${project_dir}/affinity_test"
    fi
    
    [[ -d ${target_dir} ]] && rm -r "${target_dir}"
    mkdir "${target_dir}"
    remote_target_dir=${target_dir/"${local_root}"/"${remote_root}"} # substitution
    ssh ${server_info} << EOF > '.temp'
        eval "[[ -d ${remote_target_dir} ]] && rm -r ${remote_target_dir}"
        eval "mkdir ${remote_target_dir}"
EOF
    
    # make case P16
    local remote_case_dir
    local case_dir
    local number_of_nodes=(1 1 1 1 1 1 1 2 1)
    local number_of_cores=(2 4 4 8 16 32 64 64 64)
    local bind_to_cores=(0 0 1 0 0 0 0 0 0)
    local bind_to_threads=(0 0 0 0 0 0 0 0 1)
    # todo
    local _i=0
    while ((_i<${#number_of_nodes[@]})); do
        _n=${number_of_cores[_i]}
        _N=${number_of_nodes[_i]}
        _bc=${bind_to_cores[_i]}
        _bt=${bind_to_threads[_i]}

        # deal with local files 
        case_dir="${target_dir}/N${_N}n${_n}"
        ((${_bc}==1)) && case_dir="${case_dir}bc"
        ((${_bt}==1)) && case_dir="${case_dir}bt"
        [[ -d ${case_dir} ]] && rm -r "${case_dir}"  # remove previous results
        mkdir "${case_dir}"
        eval "cp ${source_dir}/* ${case_dir}"

        # scp to remote
        remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"} # substitution
        remote_case_prm="${remote_case_dir}/case.prm"
        remote_out_dir="${remote_case_dir}/output"
        local remote_target=$(dirname "${remote_case_dir}")
        eval "${RSYNC} -r ${case_dir} ${server_info}:${remote_target}/"
    
        local status_
        while [[ true ]]; do
            ssh ${server_info} << EOF > '.temp'
                eval "[[ -e ${remote_case_prm} ]] && echo \"0\" || echo \"1\""
EOF
            status_=$(cat '.temp'| sed -n '$'p)
            ((status_==0)) && break || { cecho ${WARN} "Files haven't arrived yet, sleep for 2s"; sleep 2s; }
        done
    
        # generate job.sh file
        local addition=""
        local flag="--hold"
        ((${_bc}==1)) && flag="${flag} --bind_to=\"cores\""
        ((${_bt}==1)) && flag="${flag} --bind_to=\"threads\""
        ssh ${server_info} << EOF > ".temp"
            eval "[[ -d ${remote_out_dir} ]] || mkdir ${remote_out_dir} "
            eval "slurm.sh -N ${_N} -n ${_n} -t 24 ${addition} -P ${project} ${flag} ${remote_case_prm}"
EOF
        ((_i++))
    done
}


main(){
    ###
    # main function
    ###
    if [[ "$1" = "peloton" ]]; then
        #  affinity test on server
        #  example command line 
        #       ./aspect_lib.sh
        # todo
        set_server_info "$2"

        do_affinity_test "${project}"
    fi
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
