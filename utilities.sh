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
 
set +a  # return to default setting
