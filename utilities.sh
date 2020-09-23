#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
test_dir="${dir}/.test"  # do test in this directory
if ! [[ -d ${test_dir} ]]; then
    mkdir "${test_dir}"
fi
test_fixtures_dir="tests/integration/fixtures"

set -a  # to export every variables that will be set

################################################################################
# Unit functions
element_in() {
	# determine if element is in an array
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}


take_record() {
    # Taking record of commands
    message=$(eval echo $1)
    record_file=$(eval echo $2)

    # add time stamp to record
    timestamp () {
                echo "$(date +"%Y-%m-%d_%H-%M-%S")"

    }

    echo "$(timestamp): ${message}" >> ${record_file}
}


fix_route() {
    # substitute ~ with ${HOME}
    local dir=$(pwd)
    if [[ "${1}" =~ ^\/ ]]; then
        local output="${1}"
    elif [[ "${1}" =~ ^~\/ ]]; then
        local output=${1/#~/$HOME}
    elif [[ "${1}" =~ ^(\.\/) || "${1}" =~ ^(\.$) ]]; then
        local output=${1/#\./$dir}
    else
        local output="${dir}/${1}"
    fi
    echo "${output}"
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
get_remote_environment(){
    # get the value of a remote environment variable
    # Inputs:
    #   1: server_info
    #   2: variable name
    [[ $# = 2 ]] || return 1
    unset return_value
    local server_info=$1
    local name=$2
	ssh "${server_info}" << EOF > "${dir}/.log"
        eval "echo \${${name}}"
EOF
	return_value=$(tail -n 1 "${dir}/.log")
}


get_job_info(){
    # get info from the squeue command
    # inputs:
    #	1: job_id
    #	2: key
    unset return_value
    local _outputs
    local _temp
    if ! [[ "$1" =~ ^[0-9]*$ ]]; then
	    return 1
    fi
    _outputs=$(eval "${SQUEUE} -j ${1} 2>&1")  # thus it could be accessed through ssh
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
		[[ "${return_value}" =~ ^( )+$ ]] && return_value='NA'  # to be tested
# inorder to fix the status of the job when it just finished
	        return 0
	    fi
        ((i++))
    done
    return 2  # if the key is not find, return an error message
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
# read keys and values from a file
# Inputs:
#   $1: input file, contents for this file is:
#               key0    value0
#               key1    value1
#               ...
# Outputs:
#   keys: keys in that file
#   values: values in that file
read_keys_values(){
	local filein=$1
    # check file exist
    # check_file_exist "${filein}"
    unset keys
    unset values
	local line
	local foo

    # read in keys and values from file
	while IFS= read line; do
        IFS=' ' read -r -a foo<<< "${line}"  # construct an array from line
	    [[ -z ${keys} ]] && keys=(${foo[0]}) || keys+=(${foo[0]})
	    [[ -z ${values} ]] && values=(${foo[1]}) || values+=(${foo[1]})
    done < "${filein}"
}


read_log(){
	# read a log file
	# Inputs:
	#	$1: log file name
	local log_file=$1
	local i=0
	unset return_value0
	unset return_value1
	unset return_value2
	local line
	local foo
	while IFS= read -r line; do
        IFS=' ' read -r -a foo<<< "${line}"  # construct an array from line
        # i = 0 is the header line, ignore that
		if [[ $i -eq 1 ]]; then
			return_value0="${foo[0]}"
			return_value1="${foo[1]}"
			return_value2="${foo[2]}"
		elif [[ $i -gt 1 ]]; then
			return_value0="${return_value0} ${foo[0]}"
			return_value1="${return_value1} ${foo[1]}"
			return_value2="${return_value2} ${foo[2]}"
		fi
        ((i++))
	done < "${log_file}"
	return 0
}


write_log_header(){
	# write a header to a log file
	# Inputs:
	#	$1: log file name
	local log_file=$1
	echo "job_dir job_id ST last_time_step last_time" > "${log_file}"
}


write_log(){
    # write to a log file
    # Inputs:
    #   $1: job directory
    #   $2: job id
    #   $3: log file
    local job_dir=$1
    local job_id=$2
    local log_file=$3
    local _file
    get_job_info ${job_id} 'ST'
    quit_if_fail "get_job_info: invalid id number ${job_id} or no such stat 'ST'" 
    local ST=${return_value}

    # find the stdout file  
    for _file in ${job_dir}/*
    do
        # look for stdout file
        if [[ "${_file}" =~ ${job_id}.stdout ]]; then
            break
	fi
    done
    
    # parse stdout file
    parse_stdout ${_file}  # parse this file
    echo "${job_dir} ${job_id} ${ST} ${last_time_step} ${last_time}" >> "${log_file}"
}


clean_log(){
    # remove the record of "${case_dir}" from record
    # Inputs:
    #   $1: job directory
    #   $2: log file name
    local job_dir=$1
    local log_file=$2
    read_log "${log_file}"  # first read file
    local _find=0
    local i=2  # i=2 because we are starting from second line of a log file
    local flag=''
    for foo in ${return_value0[@]}; do
	    if [[ "${foo}" = "${job_dir}" ]]; then
		    [[ ${_find} -eq 0 ]] && { flag="${i}"; _find=1; } || flag="${flag},${i}"
	    fi
	    ((i++))
    done
    [[ ${_find} -eq 1 ]] && eval "sed -in '${flag}'d ${log_file}"  # eliminate the line of "${case_dir}"
}


################################################################################
# write time and machine time output to a file
# Inputs:
    #   $1: job directory
    #   $2: job id
    #   $3: log file
write_time_log(){
    local job_dir=$1
    local job_id=$2
    local log_file=$3

    # get machine time 
    get_job_info ${job_id} 'TIME'
    quit_if_fail "get_job_info: invalid id number ${job_id} or no such stat 'TIME'" 
    local TIME=${return_value}
    [[ -n ${TIME} && ! ${TIME} = 'NA' ]] || { cehco ${BAD} "${FUNCNAME[0]}: cannot get valid TIME for job ${job_id}"; exit 1; }

    # get CPU
    get_job_info ${job_id} 'CPU'
    quit_if_fail "get_job_info: invalid id number ${job_id} or no such stat 'CPU'" 
    local CPU=${return_value}
    [[ -n ${CPU} && ! ${CPU} = 'NA' ]] || { cehco ${BAD} "${FUNCNAME[0]}: cannot get valid CPU for job ${job_id}"; exit 1; }

    # find the stdout file  
    for _file in ${job_dir}/*
    do
        # look for stdout file
        if [[ "${_file}" =~ ${job_id}.stdout ]]; then
            break
	fi
    done
    
    # parse stdout file
    parse_stdout ${_file}  # parse this file

    # output header if file doesn't exist
    if ! [[ -e ${log_file} ]]; then 
        printf "# 1: Time step number\n" >> ${log_file}
        printf "# 2: Time\n" >> ${log_file}
        printf "# 3: Machine time\n" >> ${log_file}
        printf "# 4: CPU number\n" >> ${log_file}
    fi
    
    # output to file 
    printf "%-10s %-15s %-10s %s\n" ${last_time_step} ${last_time} ${TIME} ${CPU} >> ${log_file}

    return 0
}


################################################################################
# generate message for a tests
# Inputs:
#   $1: name of function
#   $2: number of passed tests
#   $3: number of failed tests
final_message(){
    if [[ "$3" -eq 0 ]]; then
        cecho ${GOOD} "$1: $2 tests passed"
    else
        cecho ${BAD} "$1: $2 tests passed, $3 tests failed"
    fi
}


################################################################################
# check a variable exists
# Inputs:
#   $1: name of variable
# Operations:
#   exit if variable doesn't exist
check_variable(){
    # future, should be one line
    # check variable
    echo '0'
}


################################################################################
# compare the contents of two outputs
# Inputs:
#   $1: function
#   $2: standard output
#   $3: output
compare_outputs(){
    if ! [[ $2 = "$3" ]]; then
        cecho ${BAD} "$1 failed: output - ${3} is different from standard one - ${2}"
        return 1
    fi
    return 0
}


################################################################################
# compare the contents of two files
# Inputs:
#   $1: function
#   $2: standard file
#   $3: output file
# Operations:
#   exit if file doesn't exist
compare_files(){
    # check file existence
    [[ -e $2 ]] || { cecho ${BAD} "${FUNCNAME[0]}: file $2 doesn't exist"; exit 1; }
    [[ -e $3 ]] || { cecho ${BAD} "${FUNCNAME[0]}: file $3 doesn't exist"; exit 1; }
    
    # check_file
    difference=$(diff $2 $3)
    if [[ -n ${difference} ]]; then
        cecho ${BAD} "$1 failed: output file - ${3} is different from standard one - ${2}"
        return 1
    fi
    return 0
}


################################################################################
# translate from a bash array to a python array-like output
# Inputs:
#   bash_array
# Output
#   python_array_like
bash_to_python_array()
{
    # fisrt append a '['
    python_array_like="["

    # append members in 'bash_array'
    i=0
    for member in ${bash_array[@]}; do
        # append a ',' before all but the first one
        (( i > 0 )) && python_array_like="${python_array_like}, "
        python_array_like="${python_array_like}${member}"
        ((i++))
    done

    # at last, append a ']'
    python_array_like="${python_array_like}]"
    return 0
}


################################################################################
# translate from a a python array-like string to a bash array
# currently no ' ', '[' and ']'is allowed to exist in values
# Inputs:
#   python_array_like
# Output:
#   bash_array
python_to_bash_array()
{
    # unset output
    unset bash_array

    # string substitution
    local altered_string=${python_array_like//[\[\]\ ]/}
    IFS="," read -r -a bash_array <<< "${altered_string}"
}


################################################################################
# todo
# read options from a json file
# Inputs:
#   filein: a json file
#   keys: an array of keys in json file
#   is_array: whether we want to parse an array
# Outputs:
#   value: value recorded in json file
read_json_file()
{
    # unset output variables
    unset value

    # check inputs
    [[ -e ${filein} && ${filein} =~ ".json" ]] || { cecho ${BAD} "${FUNCNAME[0]} filein must be a json file"; exit 1; }
    [[ -z ${keys} ]] && { cecho ${BAD} "${FUNCNAME[0} keys cannot be vacant"; exit 1; }
    [[ -z ${is_array} ]] && is_array="false"
    [[ ${is_array} = "true" || ${is_array} = "false" ]] || { cecho ${BAD} "${FUNCNAME[0]}: is_araay must be true or false"; exit 1; }

    # set IFS, so that we split string by " "
    IFS=" "

    # read from file
    # sed is used to get rid of "
    eval "cat ${filein} | JSON.sh | sed 's/\"//g' > ./temp"

    # construct pattern
    local i=0
    local pattern=""
    for key in ${keys[@]}; do
        (( i > 0 )) && pattern="${pattern},"
        pattern="${pattern}${key}"
        ((i++))
    done

    # grep for pattern
    local second_part=$(grep -F "[${pattern}]" temp | sed -E "s/\[${pattern}\](\t|\ )*//g")
    if [[ ${is_array} = "true" ]]; then
        python_array_like=${second_part}
        python_to_bash_array
        value=("${bash_array[@]}")
    else
        # find pattern, then split the line
        value=${second_part}
    fi

    # unset input variables
    unset is_array

    return 0
}

################################################################################
# Test functions
test_element_in(){
	local _test_array=('a' 'b' 'c d')
	if ! element_in 'a' "${_test_array[@]}"; then
		cecho ${BAD} "test_element_in failed, 'a' is not in ${_test_array}[@]"
	fi
	if element_in 'c' "${_test_array[@]}"; then
		cecho ${BAD} "test_element_in failed, 'c' is in ${_test_array}[@]"
	fi
	cecho ${GOOD} "test_element_in passed"

}


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


test_read_log(){
	local log_file="${test_fixtures_dir}/test.log"
	read_log "${log_file}"
	if ! [[ "${return_value0}" = "tests/integration/fixtures tests/integration/fixtures" && "${return_value1}" = "2009375 2009376" ]]; then
		cecho ${BAD} "test_read_log failed, return values are not correct"
		return 1
	fi
	cecho ${GOOD} "test_read_log passed"
}


test_write_log(){
    local _ofile="${test_dir}/test.log"
    if [[ -e ${_ofile} ]]; then
        # remove older file
        eval "rm ${_ofile}"
    fi
    # test 1, write a non-existent job, it should return a NA status
    write_log_header "${_ofile}"
    write_log "${test_fixtures_dir}" "2009375" "${_ofile}"
    if ! [[ -e "${_ofile}" ]]; then
        cecho ${BAD} "test_write_log fails for test1, \"${_ofile}\"  doesn't exist"
	exit 1
    fi
    _output=$(cat "${_ofile}" | sed -n '2'p)
    if ! [[ ${_output} = "${test_fixtures_dir} 2009375 NA 10 101705years" ]]
    then
        cecho ${BAD} "test_write_log fails for test2, output format is wrong"
	exit 1
    fi
    cecho ${GOOD} "test_write_log passed"

}

test_clean_log(){
	local log_source_file="${test_fixtures_dir}/test.log"
	# copy source file to log file
	local log_file="${test_dir}/test.log"
	eval "cp ${log_source_file} ${log_file}"
	# call function
	clean_log "tests/integration/fixtures" "${log_file}"
	contents=$(cat "${log_file}")  # debug
	# future compare file content
	if ! [[ "${contents}" = "job_dir job_id ST last_time_step last_time" ]]; then
		cecho ${BAD} "test_clean_log failed, file contents are not correct"
		return 1
	fi
	cecho ${GOOD} "test_clean_log passed"
}

test_fix_route() {
    local dir=$(pwd)
    # test1 test for relacing '~'
    fixed_route=$(fix_route "~/foo/ffoooo")
    [[ "${fixed_route}" = "${HOME}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 1"; exit 1; }
    # test2, test for replacing '.'
    fixed_route=$(fix_route "./foo/ffoooo")
    [[ "${fixed_route}" = "${dir}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 2"; exit 1; }
    # test3, test for replacing relative route
    fixed_route=$(fix_route "foo/ffoooo")
    [[ "${fixed_route}" = "${dir}/foo/ffoooo" ]] || { cecho ${BAD} "test_fix_route failed for test 3"; exit 1; }
    # test4, test for replacing '.'
    fixed_route=$(fix_route ".")
    [[ "${fixed_route}" = "${dir}" ]] || { cecho ${BAD} "test_fix_route failed for test 4"; exit 1; }
    # test5, test for relative address starts with .'
    fixed_route=$(fix_route ".test")
    [[ "${fixed_route}" = "${dir}/.test" ]] || { cecho ${BAD} "test_fix_route failed for test 5"; exit 1; }
    cecho ${GOOD} "test_fix_route passed"
}


main(){
	if [[ "$1" = "test" ]]; then
		# run tests by ./utilities.sh test
		test_parse_stdout
        	test_fix_route
		test_element_in
        # these two must be down with a slurm systems, future: fix it
		test_read_log
		test_write_log
		test_clean_log
	fi
}


set +a  # return to default setting


if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
