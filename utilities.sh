#!/bin/bash


set -a  # to export every variables that will be set

################################################################################
# Taking record of commands
#
# 1) Example Command:

take_record() {
    message=$(eval echo $1)
    record_file=$(eval echo $2)

    # add time stamp to record
    timestamp () {
                echo "$(date +"%Y-%m-%d_%H-%M-%S")"

    }

    echo "$(timestamp): ${message}" >> ${record_file}

}

################################################################################
# Colours for progress and error reporting
BAD="\033[1;31m"
GOOD="\033[1;32m"
WARN="\033[1;35m"
INFO="\033[1;34m"
BOLD="\033[1m"


################################################################################
# Define candi helper functions

prettify_dir() {
   # Make a directory name more readable by replacing homedir with "~"
   echo ${1/#$HOME\//~\/}
}

cecho() {
    # Display messages in a specified colour
    COL=$1; shift
    echo -e "${COL}$@\033[0m"
}

cls() {
    # clear screen
    COL=$1; shift
    echo -e "${COL}$@\033c"
}

default () {
    # Export a variable, if it is not already set
    VAR="${1%%=*}"
    VALUE="${1#*=}"
    eval "[[ \$$VAR ]] || export $VAR='$VALUE'"
}

quit_if_fail() {
    # Exit with some useful information if something goes wrong
    STATUS=$?
    if [ ${STATUS} -ne 0 ]; then
        cecho ${BAD} 'Failure with exit status:' ${STATUS}
        cecho ${BAD} 'Exit message:' $1
        exit ${STATUS}
    fi
}


################################################################################
# Define functions to work with slurm system
get_job_info(){
    # get info from the squeue command
    # inputs:
    #	1: job_id
    #	2: key
    unset return_value
    local _outputs
    local _temp
    _outputs=$(eval "squeue -j ${1} 2>&1")
    if [[ ${_outputs} =~ "slurm_load_jobs error: Invalid job id specified" ]]; then
        # catch non-exitent job id
        return_value='NA'
	return 0
    fi
    _temp=$(echo "${_outputs}" | sed -n '1'p)
    IFS=' ' read -r -a  _headers<<< "${_temp}"
    _temp=$(echo "${_outputs}" | sed -n '2'p)
    IFS=' ' read -r -a  _infos<<< "${_temp}"
    local i=0
    for element in ${_headers[@]}; do
	if [[ "$element" = "$2" ]]; then
	    return_value="${_infos[i]}"
	    return 0
	fi
	i=$i+1
    done
    return 1  # if the key is not find, return an error message
}

parse_stdout(){
	# parse from a stdout file
	# Ouputs:
	#	last_time_step(str)
	#	last_time(str): time of last time step
	local _ifile=$1
	unset last_time_step
	unset last_time
	while IFS= read -r line; do
		if [[ ${line} =~ \*\*\* ]]; then
			break
		fi
	done <<< "$(sed '1!G;h;$!d' ${_ifile})"
	last_time_step=${line#*Timestep\ }
	last_time_step=${last_time_step%:*}
	last_time=${line#*t=}
	last_time=${last_time/ /}
}

################################################################################
# Test functions
test_parse_stdout(){
	# test the parse_stdout function, return values are last timestpe and time
	local _ifile="tests/integration/fixtures/task-2009375.stdout"
	if ! [[ -e ${_ifile} ]]; then
		cecho ${BAD} "test_parse_stdout failed, no input file ${_ifile}"
		exit 1
	fi
	parse_stdout ${_ifile}  # parse this file
	if ! [[ ${last_time_step} = "10" ]]; then
		cecho ${BAD} "test_parse_stdout failed, time_step is wrong"
		exit 1
	fi
	if ! [[ ${last_time} = "101705years" ]]; then
		cecho ${BAD} "test_parse_stdout failed, time is wrong"
		exit 1
	fi
	cecho ${GOOD} "test_parse_stdout passed"
}

if [[ "$1" = "test" ]]; then
	# run tests by ./utilities.sh test
	test_parse_stdout
fi
 
set +a  # return to default setting
