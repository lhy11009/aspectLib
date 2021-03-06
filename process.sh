#!/bin/bash
# todo_future : put tests to seperate files and have a test_manager

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source "${dir}/utilities.sh"

test_dir="${dir}/.test"  # do test in this directory
if ! [[ -d ${test_dir} ]]; then
    mkdir "${test_dir}"
fi
test_fixtures_dir="${dir}/tests/integration/fixtures"


usage(){
	# future
	echo "foo"
}

################################################################################
# parse_options
parse_options(){
    # parse parameters from command line
    # future pass in name of case
    while [ -n "$1" ]; do
      param="$1"
      case $param in
        -h|--help)
          echo ""  # help information
          exit 0
        ;;
        #####################################
        # filename
        #####################################
        [^-]*)
          filename=$(fix_route "$param")
        ;;
        #####################################
        # number of total tasks
        #####################################
        -n)
          shift
          total_tasks="${1}"
        ;;
        -n=*|--total_tasks=*)
          total_tasks="${param#*=}"
        ;;
        #####################################
        # number of nodes
        #####################################
        -N)
          shift
          nnode="${1}"
        ;;
        -N=*|--nnode=*)
          nnode="${param#*=}"
        ;;
        #####################################
        # time in hour
        #####################################
        -t)
          shift
          time_by_hour="${1}"
        ;;
        -t=*|--time=*)
          time_by_hour="${param#*=}"
        ;;
        #####################################
        # partition
        #####################################
        -p)
          shift
          partition="${1}"
        ;;
        -p=*|--partition=*)
          partition="${param#*=}"
        ;;
        #####################################
        # memory per cpu
        #####################################
        -m)
          shift
          mem_per_cpu="${1}"
        ;;
        -m=*|--mem-per-cup=*)
          mem_per_cpu="${param#*=}"
        ;;
        #####################################
        # log file
        #####################################
        -l)
          shift
          local temp="${1}"
	  log_file=$(fix_route "${temp}")
        ;;
        -l=*|--log_file=*)
          local temp="${param#*=}"
	  log_file=$(fix_route "${temp}")
        ;;
        #####################################
        # log file for time
        #####################################
        -lt)
          shift
          local temp="${1}"
	  log_file_time=$(fix_route "${temp}")
        ;;
        -lt=*|--log_file_time=*)
          local temp="${param#*=}"
	  log_file_time=$(fix_route "${temp}")
        ;;
        #####################################
        # project
        #####################################
        -P)
          shift
          project="${1}"
        ;;
        -P=*|--project=*)
          project="${param#*=}"
        ;;
        #####################################
        # bool value
        #####################################
        -b)
          shift
          bool="${1}"
        ;;
        -b=*|--bool=*)
          bool="${param#*=}"
        ;;
        #####################################
        # float value
        #####################################
        -f)
          shift
          float="${1}"
        ;;
        -f=*|--float=*)
          float="${param#*=}"
        ;;
      esac
      shift
    done
    
    # check values
    [[ -z ${bool} || ${bool} = "true" || ${bool} = "false" ]] || { cecho ${BAD} "${FUNCNAME[0]}: bool value must be true or false"; exit 1; }
    [[ -z ${float} || ${float} =~ ^[0-9\.]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: entry for \${float} must be a float value"; exit 1; }
}

################################################################################
# submit case
submit(){
	_sbatch=$(which sbatch)
	if [[ -z ${_sbatch} ]]; then
		_mpirun=$(which mpirun)
		[[ -z ${_mpirun} ]] && { cecho $BAD "${FUNCNAME[0]}: this is neither slurm or mpi in system"; exit 1; } || submit_mpirun
	else
		submit_slurm
	fi
}

################################################################################
# submit case through slurm
submit_slurm(){
	# todo
	echo "foo"
}

################################################################################
# submit case through mpirun
# Inputs:
# 	$1: name of project(optional)
submit_mpirun(){

	# todo
	project="$1"

	# fix executable
	[[ -z ${project} ]] && executable="${ASPECT_SOURCE_DIR}/build/aspect" || executable="${ASPECT_SOURCE_DIR}/build_${project}/aspect"

	_dir=$(dirname "${filename}"); 
	_base=$(basename "${filename}"); 
	cd "${_dir}"
	mpirun -np "$total_tasks" "$executable" "${basename}" >"job.stdout" 2>"job.stderr" &
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
		((i++))
	done
	# show content
	echo "content of log file:"
	eval "cat ${log_file}"
}

clean_NA(){
	local log_file=$1
	local job_dirs
	local job_ids
	local job_dir
	local job_id
	local job_ST
	# look into log file
	read_log "${log_file}"
	# write log file
	write_log_header "${log_file}"
	IFS=' ' read -r -a job_dirs <<< ${return_value0}
	IFS=' ' read -r -a job_ids <<< ${return_value1}
	IFS=' ' read -r -a job_STs <<< ${return_value2}
	i=0
	for job_id in ${job_ids[@]}; do
		job_dir="${job_dirs[i]}"
		job_ST="${job_STs[i]}"
		[[ "${job_ST}" = 'NA' ]] || write_log "${job_dir}" "${job_id}" "${log_file}"
		((i++))
	done
	# updata log file
}


update_from_server(){
	# update log file from remote
	local server_info=$1
	local remote_log_file=$2
	local local_log_file=$3
    # use ssh to update log file on server side, \$ escapes '$' so that it is called on the remote side
    ssh "$server_info" eval "\${ASPECT_LAB_DIR}/process.sh update ${remote_log_file}"
    # use rsync to update log file on local side
    eval "${RSYNC} ${server_info}:${remote_log_file} ${local_log_file}"
}

clean_NA_from_server(){
	# update log file from remote
	local server_info=$1
	local remote_log_file=$2
	local local_log_file=$3
    # use ssh to update log file on server side, \$ escapes '$' so that it is called on the remote side
    ssh "$server_info" eval "\${ASPECT_LAB_DIR}/process.sh clean ${remote_log_file}"
    # use rsync to update log file on local side
    eval "${RSYNC} ${server_info}:${remote_log_file} ${local_log_file}"
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
		echo "${RSYNC} -avur ${server_info}:${job_dir}/* ${target_dir}/"
		eval "${RSYNC} -avur ${server_info}:${job_dir}/* ${target_dir}/"
	done
}

################################################################################
# show local jobs ran with mpirun
show_local_jobs(){
   	eval "ps -aux | grep build.*aspect" # show all process
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

test_clean_NA(){
	local log_file="${test_dir}/test.log"
	local correct_output_file="${test_fixtures_dir}/outputs/process_sh_test_clean_NA.output"
	eval "cp ${test_fixtures_dir}/test.log ${log_file}"
	clean_NA "${log_file}"
	if [[ $? -eq 1 ]]; then
		cecho ${BAD} "test_update failed, error reading log file ${log_file}"
		exit 1
	fi
	if ! [[ $(diff "${log_file}" "${correct_output_file}") = '' ]]; then
		cecho ${BAD} "test_update failed, output format is wrong"
		exit 1
	fi
	cecho ${GOOD} "test_clean_NA passed"
}

test_update_from_server(){
    # todo_future: have a function to check all the environmental variables on server
    # test 1: call update_from_server function
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
	remote_test_log_file="${return_value}/.test/test.log"
    # mannually fix the test through a ssh session
    # the sed line is tricky, it is used to mannually fit the route in the test
    # otherwise it could be mischievous on different machines
    ssh "${server_info}" << EOF
        [[ -e "${remote_test_log_file}" ]] && eval "rm ${remote_test_log_file}"
        eval "cp ${remote_log_file} ${remote_test_log_file}"
        eval "sed -i \"2,3s+^+${return_value}\\/+g\" \"${remote_test_log_file}\""
EOF
	update_from_server "${server_info}" "${remote_test_log_file}" "${local_log_file}"
    if ! [[ -e "${local_log_file}" ]]; then
		cecho ${BAD} "test_update failed, output file doesn't exist"
		exit 1
    fi
    local line2=$(sed -n '2'p "${local_log_file}")  # second line of local log file
    local line3=$(sed -n '3'p "${local_log_file}")  # third line of local log file
    local correct_line2="${return_value}/tests/integration/fixtures 2009375 NA 10 101705years"
    local correct_line3="${return_value}/tests/integration/fixtures 2009376 NA 11 101706years"
    if ! [[ "${line2}" = "${correct_line2}" && "${line3}" = "${correct_line3}" ]]; then
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
	_command="$1"
    parse_options "$@"  # todo, parse option with '-'

	if [[ "${_command}" = "test" ]]; then
		test_update
		test_clean_NA

	elif [[ "${_command}" = "remote_test" ]]; then
        # test run script on remote
		if ! [[ $# -eq 2 ]]; then
			cecho ${BAD} "with \"remote_test\" command, \$2 must be given for server_info"
            exit 1
		fi
        set_server_info "$2"
		test_update_from_server "${server_info}"
		test_update_outputs_from_server "${server_info}"

	elif [[ "${_command}" = "add" ]]; then
		if [[ "$#" -eq 4 ]]; then
			cecho ${BAD} "with \"update\" command, \$2 \$3 and \$4 must be given for job_dir, job_id, log_file"
            		exit 1
	    	fi
		local job_dir="$2"
		local job_id="$3"
		local log_file=$4
		write_log "${job_dir}" "${job_id}" "${log_file}"

	elif [[ "${_command}" = "update" ]]; then
        # todo_future, strip root dir from output dir
		local log_file=$2
		if [[ "${log_file}" = '' ]]; then
			cecho ${BAD} "with \"update\" command, a \$2 must be given for the log_file variable"
            exit 1
		fi
	       	update "${log_file}"

	elif [[ "${_command}" = "clean" ]]; then
        # todo_future, strip root dir from output dir
		local log_file=$2
		if [[ "${log_file}" = '' ]]; then
			cecho ${BAD} "with \"clean\" command, a \$2 must be given for the log_file variable"
            exit 1
		fi
	       	clean_NA "${log_file}"

	elif [[ "${_command}" = "update_from_server" ]]; then
        # download new log file from server
        # todo_future, use a config file for configration
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"update_from_server\" command, \$2, \$3 must be given for server information, local log file"
            exit 1
		fi
		local local_log_file=$2
		set_server_info "$3"
        # fix route
        local_log_file=$(fix_route "${local_log_file}")
        # figure out remote directory
        get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	    remote_log_file=${local_log_file/"${dir}"/"${return_value}"}
        # call function to transfer file
		update_from_server "${server_info}" "${local_log_file}" "${remote_log_file}"

	elif [[ "${_command}" = "clean_from_server" ]]; then
        # download new log file from server
        # todo_future, use a config file for configration
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"clean_from_server\" command, \$2, \$3 must be given for server information, local log file"
            exit 1
		fi
		local local_log_file=$2
		set_server_info "$3"
        # fix route
        local_log_file=$(fix_route "${local_log_file}")
        # figure out remote directory
        get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	    remote_log_file=${local_log_file/"${dir}"/"${return_value}"}
        # call function to transfer file
	    clean_NA_from_server "${server_info}" "${local_log_file}" "${remote_log_file}"

	elif [[ "${_command}" = "update_outputs_from_server" ]]; then
        # test update_outputs_from_server
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"update_outputs_from_server\" command, \$2 and \$3 must be given for server_info and log_files on local side"
            exit 1
		fi
		local local_log_file=$2
		set_server_info "$3"
        update_outputs_from_server "${server_info}" "${local_log_file}"

	elif [[ "${_command}" = "remove" ]]; then
        # remove both local and remote log file
		if ! [[ $# -eq 3 ]]; then
			cecho ${BAD} "with \"clean\" command, \$2, \$3 must be given for local log file and server information"
            exit 1
		fi
		local local_file=$2
		set_server_info "$3"
        # fix route
        local_file=$(fix_route "${local_file}")
        # figure out remote directory
        get_remote_environment ${server_info} "ASPECT_LAB_DIR"
	    remote_file=${local_file/"${dir}"/"${return_value}"}
        # remove local and remote files
        [[ -e "${local_file}" ]] && eval "rm  -r ${local_file}"
        ssh "${server_info}" eval "[[ -e "${remote_file}" ]] && rm -r ${remote_file}"

	elif [[ "${_command}" = "show_local_jobs" ]]; then
        # show local jobs ran with mpirun
        #   example command lines:
		#		./process.sh show_local_jobs
		show_local_jobs
    
	elif [[ ${_command} = "submit" ]]; then
		# submit jobs
		submit

	elif [[ "${_command}" = '-h' || "${_command}" = '--help' ]]; then
		usage

	else
		cecho ${BAD} "bad command \"${_command}\""
		usage
	fi
}


if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
