#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source "${dir}/utilities.sh"
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
	# look into log file
	eval "slurm.sh"  # todo
	# updata log file
	eval "slurm.sh" # todo
}


update_from_server(){
	local server_info=$2
	local local_log_file=$3
	local remote_log_file=$4
	# update log file
	echo "rsync -avur "  # todo
}


update_outputs_from_server(){
	# copy case output from server
	local server_info=$2
	local local_log_file=$3
	local line
	for IFS=read -r line; do
		echo "rsync -avur"  # todo
	done <<< $(cat "${log_file}")
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
	echo "foo" # todo
}


test_update_outputs_from_server(){
	echo "foo" # todo
}


################################################################################
# main function
################################################################################
main(){
	if [[ "$1" = "test" ]]; then
		test_update
	elif [[ "$1" = "remote_test" ]]; then
		test_update_from_server
		test_update_outputs_from_server
	elif [[ "$1" = "update" ]]; then
		local log_file=$2
		if [[ "${log_file}" = '' ]]; then
			cecho ${BAD} "with \"update\" command, a $2 must be given for the log_file variable"
		fi
	       	update "${log_file}"
	elif [[ "$1" = "update_from_server" ]]; then
		local server_info=$2
		local local_log_file=$3
		local remote_log_file=$4
		if [[ "${local_log_file}" = '' || "${remote_log_file}" = '' ]]; then
			cecho ${BAD} "with \"update_from_server\" command, $3and $4 must be given for log_files on local side and log_files on remote side"
		fi
		update_from_server "${server_info}" "${local_log_file}" "${remote_log_file}"
	elif [[ "$1" = "update_outputs_from_server" ]]; then
		local server_info=$2
		local local_log_file=$3
		if [[ "${remote_log_file}" = '' ]]; then
			cecho ${BAD} "with \"update_from_server\" command, $3 must be given for log_files on local side"
		fi
		update_outputs_from_server "${server_info}" "${remote_log_file}"
	elif [[ "$1" = '-h' || "$1" = '--help' ]]; then
		usage
	else
		cecho ${BAD} "bad command \"$1\""
		usage
	fi
}


main
