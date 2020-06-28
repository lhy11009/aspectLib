#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source "${dir}/utilities.sh"
source "${dir}/extra_configurations.sh"

test_dir="${dir}/.test"  # do test in this directory
if ! [[ -d ${test_dir} ]]; then
    mkdir "${test_dir}"
fi
test_fixtures_dir="${dir}/tests/integration/fixtures"


usage(){
	#todo
	echo "foo"
}


update(){
	local log_file=$1
	local job_dirs
	local job_ids
	local job_dir
	local job_id
	# look into log file
	read_log "${log_file}"
	# write log file
	write_log_header "${log_file}"
	IFS=' ' read -r -a job_dirs <<< ${return_value0}
	IFS=' ' read -r -a job_ids <<< ${return_value1}
	i=0
	for job_id in ${job_ids[@]}; do
		job_dir="${job_dirs[i]}"
		write_log "${job_dir}" "${job_id}" "${log_file}"
		i=$i+1
	done
	# updata log file
}


update_from_server(){
	# update log file from remote
    [[ ${RSYNC} = '' ]] && local RSYNC='rsync'
	local server_info=$1
	local remote_log_file=$2
	local local_log_file=$3
    eval "${RSYNC} -v ${server_info}:${remote_log_file} ${local_log_file}"
}


update_outputs_from_server(){
	# copy case output from server
    # files to update is read from a local log file, so before calling this function
    # a local log file should be generated;
    # local directory and remote directory are withing the same project, so the
    # path for the project is substituted from the remote directory to the
    # local directory
    # Inputs:
    #   $1: server_info, user@address
    #   $2: local_log_file
    local project="TwoDSubduction"  # todo: loop for different project
    [[ ${RSYNC} = '' ]] && local RSYNC='rsync' # RSYNC could be set to substitute rsync,todo: look for better solution
	local server_info=$1
	local local_log_file=$2
    local local_root=$(eval "echo \${${project}_DIR}")  # get local root dir
    local target_dir
    get_remote_environment "${server_info}" "${project}_DIR"  # get remote dir
    local remote_root=${return_value}
	read_log "${local_log_file}"
    local job_dirs  # read local log file
	IFS=' ' read -r -a job_dirs <<< ${return_value0}
	for job_dir in ${job_dirs[@]}; do
        target_dir=${job_dir/"${remote_root}"/"${local_root}"}
        [[ -d ${target_dir} ]] || mkdir -p "${target_dir}"  # mkdir target dir if it doesn't exist
		eval "${RSYNC} -avur ${server_info}:${job_dir}/* ${target_dir}/"
	done
}
################################################################################
# test functions
################################################################################
test_update(){
	local log_file="${test_dir}/test.log"
	local correct_output_file="${test_fixtures_dir}/outputs/process_sh_test_update.output"
	eval "cp ${test_fixtures_dir}/test.log ${log_file}"
	update "${log_file}"
	if [[ $? -eq 1 ]]; then
		cecho ${BAD} "test_update failed, error reading log file ${log_file}"
		exit 1
	fi
	if ! [[ $(diff "${log_file}" "${correct_output_file}") = '' ]]; then
		cecho ${BAD} "test_update failed, output format is wrong"
		exit 1
	fi
	cecho ${GOOD} "test_update passed"
}


test_update_from_server(){
    # todo_future: have a function to check all the environmental variables on server
	local correct_output_file="${test_fixtures_dir}/test.log"
	local server_info=$1
	local local_log_file="${test_dir}/test.log"
	if [[ -e "${local_log_file}" ]]; then
		# remove old file
		eval "rm ${local_log_file}"
	fi
    # figure out remove file
    # todo_future: fix the return value
    get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	remote_log_file="${return_value}/tests/integration/fixtures/test.log"
	update_from_server "${server_info}" "${remote_log_file}" "${local_log_file}"
    if ! [[ -e "${local_log_file}" ]]; then
		cecho ${BAD} "test_update failed, output file doesn't exist"
		exit 1
    fi
	if ! [[ $(diff "${local_log_file}" "${correct_output_file}") = '' ]]; then
		cecho ${BAD} "test_update failed, output format is wrong"
		exit 1
	fi
	cecho ${GOOD} "test_update_from_server passed"
}


test_update_outputs_from_server(){
    # todo: this is a test need to be done in the project folder
    # so there is need to make this test when installing this bundle
    # Inputs:
    #   $1:server_info
    server_info=$1
    local project="TwoDSubduction"  # todo: loop for different project
    local local_root=$(eval "echo \${${project}_DIR}")  # get local root dir
    local target_dir="${local_root}/tests/update_outputs_from_server_tests"
    [[ -d "${target_dir}" ]] && rm -r "${target_dir}"  # remove old dir
    local_log_file="${local_root}/tests/test.log"
    update_outputs_from_server "${server_info}" "${local_log_file}"
    if ! [[ -e "${target_dir}/foo1/foo1a" && -e "${target_dir}/foo1/foo1b" && -e "${target_dir}/foo1/foo1c" ]]; then
        cecho ${BAD} "test_update_outputs_from_server failed, file in foo1 doesn't exist";
        return 1
    fi
    if ! [[ -e "${target_dir}/foo2/foo2a" && -e "${target_dir}/foo2/foo2b" && -e "${target_dir}/foo2/foo2c" ]]; then
        cecho ${BAD} "test_update_outputs_from_server failed, file in foo2 doesn't exist";
        return 1
    fi
    cecho ${GOOD} "test_update_outputs_from_server succeed"
}

################################################################################
# main function
################################################################################
main(){
	if [[ "$1" = "test" ]]; then
		test_update
	elif [[ "$1" = "remote_test" ]]; then
        # test run script on remote
		if ! [[ $# -eq 2 ]]; then
			cecho ${BAD} "with \"remote_test\" command, \$2 must be given for server_info"
            exit 1
		fi
        local server_info=$2
		test_update_from_server "${server_info}"
		test_update_outputs_from_server "${server_info}"
	elif [[ "$1" = "update" ]]; then
		local log_file=$2
		if [[ "${log_file}" = '' ]]; then
			cecho ${BAD} "with \"update\" command, a \$2 must be given for the log_file variable"
            exit 1
		fi
	       	update "${log_file}"
	elif [[ "$1" = "update_from_server" ]]; then
        # download new log file from server
        # todo_future, use a config file for configration
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"update_from_server\" command, \$2, \$3 must be given for server information, local log file"
            exit 1
		fi
		local local_log_file=$2
		local server_info=$3
        # figure out remote directory
        get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	    remote_log_file=${local_log_file/"${dir}"/"${return_value}"}
		update_from_server "${server_info}" "${local_log_file}" "${remote_log_file}"
	elif [[ "$1" = "update_outputs_from_server" ]]; then
        # test update_outputs_from_server
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"update_outputs_from_server\" command, \$2 and \$3 must be given for server_info and log_files on local side"
            exit 1
		fi
		local server_info=$2
		local local_log_file=$3
		update_outputs_from_server "${server_info}" "${remote_log_file}"
	elif [[ "$1" = '-h' || "$1" = '--help' ]]; then
		usage
	else
		cecho ${BAD} "bad command \"$1\""
		usage
	fi
}


if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
