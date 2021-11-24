#!/bin/bash

################################################################################
# Parse block outputs from a log file that ASPECT generates
#
# Dependencies:
#
# Example Usage:
#    analyse affinity results:
#        ./bash_scripts/parse_block_output.sh analyze_affinity_test_results 
# /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/results/spherical_shell_expensive_solver/peloton-ii-32tasks-core-openmpi-4.0.1/output_16_2_1 
#    
#    case runtime info:
#        ./parse_block_output.sh case_runtime /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/bsq
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${ASPECT_LAB_DIR}/utilities/bash_scripts/utilities.sh"


usage(){
  # usage of this script
    _text="
${BASH_SOURCE[0]}

Parse job run time information

Dependencies:
   env:
       ASPECT_LAB_DIR
       ASPECT_SOURCE_DIR

Example Usage:
   case run time info (last step):
       
       Lib_parse_case.sh case_info /group/billengrp-mpi-io/lochy/TwoDSubduction/wb/wb_cart_4
   
   case run time info by given an id (must be running):

      Lib_parse_case.sh case_info_with_id 250256 
"
    printf "${_text}"

}


parse_block_output(){
    ##
    # Parse value in output, looking for the block output of aspect
    # Inputs:
    #   $1(str): logfile
    # Outputs:
    #   ??: 
    #       entries:
    #           0: Total wallclock time elapsed since start
    #           1: Assemble Stokes system
    #           2: Assemble composition system
    #           3: Assemble temperature system
    #           4: Build Stokes preconditioner
    #           5: Build composition preconditioner
    #           6: Build temperature preconditioner
    #           7: Initialization
    #           8: Postprocessing
    #           9: Setup dof systems
    #           10: Setup initial conditions
    #           11: Setup matrices
    #           12: Solve Stokes system
    #           13: Solve composition system
    #           14: Solve temperature system
    unset return_values
    local logfile="$1"
    local key="$2"
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile doesn't exist"
    # read file from the end
    local parse_results=$(eval "awk '/${key}/{print}' ${logfile} | sed 's/${key}//g' |  awk '{print \$5}' | sed ':a;N;\$!ba;s/\n/ /g'")
    IFS=' '; return_values=(${parse_results})
}


parse_block_output_wallclock(){
    ##
    # Parse value in output, looking for Total wallclock time
    # It turns output the previous one doesn't work for wallclock time
    unset return_values
    local logfile="$1"
    local key="Total wallclock time elapsed since start"
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile doesn't exist"
    # read file from the end
    local parse_results=$(eval "awk '/${key}/{print}' ${logfile} | sed 's/${key}//g' |  awk '{print \$3}' | sed ':a;N;\$!ba;s/\n/ /g'")
    IFS=' '; return_values=(${parse_results})
}


parse_block_output_to_file(){
    ##
    # Parse value in output, looking for the block output of aspect
    # Inputs:
    #   $1(str): logfile
    #   $2(str): ofile
    #   $3-: keys
    local logfile="$1"
    local ofile="$2"
    # checkfile exist
    [[ -e ${logfile} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile(${logfile}) doesn't exist"

    local key="Total wallclock time elapsed since start";
    local header="# ${key}\n"
    parse_block_output_wallclock "${log_file}"
    output=$(echo ${return_values[@]} | sed -E "s/[^0-9.e+]/ /g")  # could be wrong when it is scientific expression
    local contents="${output}"


    # loop for key
    local i=3
    key=${!i}
    while [[ -n ${key} ]]; do
        header="${header}# ${key}\n"
        parse_block_output "${log_file}" "${key}"
        output=$(echo ${return_values[@]} | sed -E "s/[^0-9.e+]/ /g")  # could be wrong when it is scientific expression
        [[ -n ${contents} ]] && contents="${contents}\n${output}" || contents=${output}
        ((i++))
        key=${!i}
    done

    # output
    printf "${header}" > "${ofile}"
    printf "${contents}" >> "${ofile}"
}

################################################################################
# functions to parse run time info
################################################################################

parse_run_time_info_last_step(){
    # parse run time output
    # $1(str): log file
    [[ -e ${1} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile(${1}) doesn't exist"
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output" "$1")
    IFS=$'\n'; local entries=(${outputs})
    IFS=' '; return_value=(${entries[*]: -1})
    printf "Step\t Time(unit defined in prm)\t Wall_clock(s)\n"
    printf "${return_value[*]}\n"
}

parse_run_time_info_last_step_with_id(){
    # parse run time output with job id
    # $1(str): job id
    printf "job_id: $1\n"
    locate_workdir_with_id "$1"
    [[ $? == 1 ]] && return 1
    [[ -d "${return_value}" ]] ||  { cecho ${BAD} "${FUNCNAME[0]}: work directory(${return_value}) doesn't exist"; exit 1; }
    local log_file="${return_value}/output/log.txt"
    [[ -e "${log_file}" ]] || { cecho "${BAD}" "${FUNCNAME[0]}: log file (${log_file}) doesn't exist"; exit 1; }
    parse_run_time_info_last_step "${log_file}"
    printf "\n"
    return 0
}

parse_all_time_info(){
    # show all jobs info
    local outputs=$(squeue -u "${USER}" -o %A)
    IFS=$'\n'; local job_ids=(${outputs})
    for job_id in ${job_ids[@]}; do
        [[ ${job_id} =~ ^[0-9]+$ ]] && parse_run_time_info_last_step_with_id ${job_id}
    done
}
 
locate_workdir_with_id(){
    # locate work directory of a job with job id
    # $1 (str): job
    local outputs=$(scontrol show job "$1" | grep WorkDir)
    if [[ -n ${outputs} ]]; then
        IFS='='; local entries=(${outputs})
	return_value=${entries[*]: -1}
	printf "workdir: ${return_value}\n"
	return 0
    else
        cecho "${BAD}" "no job (id: $1) found"; return 1
    fi
}

################################################################################
# functions to submit cases
################################################################################

submit_case_peloton_rome(){
    ####
    # todo
    # submit case to slurm
    # Inputs:
    #   $1: case directory
    now="${pwd}"
    cd "$1"
    eval "sbatch -A biilen job_rome.sh"
    cd "${now}"
    return 0
}

restart_case(){
    ####
    # todo
    # restart a case
    # Inputs:
    #   $1: case directory
    ####
    # change options in the case.prm file
    local case_dir="$1"
    local prm_file="${case_dir}/case.prm"
    [[ -e ${prm_file} ]] || { cecho "${BAD}" "No such file ${prm_file}"; exit 1; } 
    util_rewrite_prm "${prm_file}" "Resume Computation" "true"
    submit_case_peloton_rome "${case_dir}"
    # todo
    return 0
}

check_time_restart_case(){
    ####
    # todo
    # restart a case if the run time of that case is not reached
    # Inputs:
    #   $1: case directory
    #   $2: run time
    ####
    local case_dir="$1"
    local run_time_std="$2"
    # read run_time
    local log_file="$1/output/log.txt"
    [[ -e ${log_fiel} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile(${log_file}) doesn't exist"
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output" "${log_file}")
    IFS=$'\n'; local entries=(${outputs})
    IFS=' '; return_value=(${entries[*]: -1})
    # run_time = 
    # [[ ]] && restart_case "$1"
    return 0
}

# future remove unneeded functions and add one to parse case info(i.e. id job given, case name, wall clock, core*hrs, step

main(){
    if [[ "$1" = "-h" ]]; then
        usage
    elif [[ "$1" = "parse_block_results" ]]; then
        # this doesn't work for $3 with whitespace in it.
	    [[ -n "$2" ]] || { cecho "${BAD}" "no log file given (\$2)"; exit 1; }
	    [[ -n "$3" ]] || { cecho "${BAD}" "no key given (\$3)"; exit 1; }
        parse_block_output "$2" "$3"
        printf "${return_values}"
    
    elif [[ "$1" = "analyze_affinity_test_results" ]]; then
	    # analyze affinity test results
        local log_file="$2"
        local ofile="$3"
        parse_block_output_to_file "${log_file}" "${ofile}" "Assemble Stokes system" "Solve Stokes system"
    elif [[ "$1" = "case_info" ]]; then
        # todo
	[[ -n "$2" ]] || { cecho "${BAD}" "no case directory given (\$2)"; exit 1; }
	[[ -d "$2" ]] || { cecho "${BAD}" "directory given (\$2) doesn't exist"; exit 1; }
	log_file="$2/output/log.txt"
	[[ -e "${log_file}" ]] || { cecho "${BAD}" "log file (${log_file}) doesn't exist"; exit 1; }
	parse_run_time_info_last_step "${log_file}"
    elif [[ "$1" = "case_info_with_id" ]]; then
	[[ -n "$2" ]] || { cecho "${BAD}" "no id number given (\$2)"; exit 1; }
	parse_run_time_info_last_step_with_id "$2"
    elif [[ "$1" = "all_case_info" ]]; then
	parse_all_time_info
    else
    	cecho "${BAD}" "option ${1} is not valid\n"
    fi
}

set +a  # return to default setting

if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
