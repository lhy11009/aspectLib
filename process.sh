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
    local project="TwoDSubduction"  # todo: loop for different project
	local server_info=$1
	local local_log_file=$2
    local local_root=$(eval "echo \${${project}_DIR}")  # get local root dir
	local line
	read_log "${local_log_file}"
	IFS=' ' read -r -a job_dirs <<< ${return_value0}
	for job_dir in ${job_dirs[@]}; do
		echo "rsync -avur "  # todo
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
	local correct_output_file="${test_fixtures_dir}/test.log"
	local server_info=$1
	local local_log_file="${test_dir}/test.log"
	if [[ -e "${local_log_file}" ]]; then
		# remove old file
		eval "rm ${local_log_file}"
	fi
    # figure out remove file
    get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	echo "REMOTE_ASPECT_LAB_DIR: ${return_value}"  # screen output
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
		if ! [[ $# -eq 4 ]]; then
			cecho ${BAD} "with \"update_from_server\" command, \$2, \$3 and \$4 must be given for server information, log_files on local side and log_files on remote side"
            exit 1
		fi
		local server_info=$2
		local local_log_file=$3
		local remote_log_file=$4
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
