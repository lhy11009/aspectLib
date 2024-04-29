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
       
       Lib_parse_case case_info /group/billengrp-mpi-io/lochy/TwoDSubduction/wb/wb_cart_4
   
   case run time info by given an id (must be running):

      Lib_parse_case case_info_with_id 250256 

   information about all the runnign case:
      
      Lib_parse_case all_case_info
    
   export run time info of a case to a log file
    
      Lib_parse_case export_case_info ~/ASPECT_PROJECT/TwoDSubduction/wb_sd_issue_2/wb_sph_cdd50_substract_T_op40_20Ma_hr export_test.txt
    
   export run time info with a directory (all the cases in this directory) to \"cases.log\":
      Lib_parse_case export_all_case_info_in_directory . 

   Restart case by checkng running time:
      Lib_parse_case check_time_restart /group/billengrp-mpi-io/lochy/TwoDSubduction/wb/wb_cart_4 10e6

   Restart cases by checking run time in a directory
      Lib_parse_case check_time_restart_directory ./EBA_CDPT 60e6 high2

   Copy the case files to a new location, including only one step in the output
      (the last number is the visualization step (step + initiation adaptive refinement level))
      Lib_parse_case copy_case_output_by_vtu_snapshot ~/ASPECT_PROJECT/aspectLib/tests/integration/big_fixtures/test_bash_parse/eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr9_yd100
      .test/test_base_parse 0
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
    IFS=$'\n'; local entries=(${outputs})  # each member in entries is a separate line
    IFS=' '; return_value=(${entries[*]: -1})  # get the last line and split with ' '
    [[ ${return_value} =~ "#" ]] && return_value=""  # if only header is retured, reset
    printf "Step\t Time(unit defined in prm)\t Wall_clock(s)\n"
    printf "${return_value[*]}\n"
}

parse_run_time_info_last_step_in_dir(){
    # parse run time output from cases in a directory
    [[ -d ${1} ]] || cecho ${BAD} "${FUNCNAME[0]}: diretory path(${1}) doesn't exist"
    local dir="$1"
    local log_file
    for sub_dir in "${dir}"/*; do
        log_file="${sub_dir}/output/log.txt"
	if [[ -e "${log_file}" ]]; then
	    local case_name=$(basename "${sub_dir}")
	    echo "${case_name}"
	    parse_run_time_info_last_step "${log_file}"
	fi
    done
}

parse_run_time_info_last_step_with_id(){
    # parse run time output with job id
    # $1(str): job id
    printf "job_id: $1\n"
    locate_workdir_with_id "$1"
    [[ $? == 1 ]] && return 1
    [[ -d "${return_value}" ]] ||  { cecho ${BAD} "${FUNCNAME[0]}: work directory(${return_value}) doesn't exist"; exit 1; }
    local log_file="${return_value}/output/log.txt"
    if  [[ -e "${log_file}" ]]; then
    	parse_run_time_info_last_step "${log_file}"
    else
	cecho "${WARNING}" "${FUNCNAME[0]}: log file (${log_file}) doesn't exist.\
It's either this case just started or the file is lost."
    fi  
    printf "\n"
    return 0
}

parse_export_run_time_info(){
    # parse run time info and export to a file
    # $1 (str): case directory
    # $2 (str): file to export
    # check
    local temp
    local case_dir="$1"
    local case_name=$(basename "${case_dir}")
    data0=("${case_name}") # case name
    local file_path="$2"
    local append="$3"  # append to old file
    local log_file="${case_dir}/output/log.txt"
    [[ -d "${case_dir}" ]] ||  { cecho ${BAD} "${FUNCNAME[0]}: case directory(${case_dir}) doesn't exist"; exit 1; }
    headers=("case_name" "status" "step" "last_restart_step" "Time" "wallclock" "total_wallclock")
    # parse case info
    [[ -e ${log_file} ]] || cecho ${BAD} "${FUNCNAME[0]}: logfile(${log_file}) doesn't exist"
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_block_output" "${log_file}")
    IFS=$'\n'; local entries=(${outputs})  # each member in entries is a separate line
    IFS=' '; return_value=(${entries[*]: -1})  # get the last line and split with ' '
    temp=$(util_neat_word "${return_value[0]}"); data2=("${temp}")  # step
    temp=$(util_neat_word "${return_value[1]}"); data4=("${temp}")  # time
    # status
    temp=$(get_case_status "${case_dir}"); data1=("${temp}")
    # parse resume-computation info
    local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_resume_computation" "${log_file}")
    IFS=$'\n'; local entries=(${outputs})  # each member in entries is a separate line
    IFS=' '; return_value=(${entries[*]: -1})  # get the last line and split with ' '
    temp=$(util_neat_word "${return_value[2]}"); data5=("${temp}")  # wall clock since last restart
    IFS=' '; return_value=(${entries[*]: -2})  # get the previous line and split with ' '
    temp=$(util_neat_word "${return_value[0]}"); data3=("${temp}")  # step of last restart
    local total='0.0'
    for entry in "${entries[@]}"; do
        return_value=(${entry})
        temp="${return_value[2]}"
        total=$(awk -v a="${total}" -v b="${temp}" 'BEGIN{printf "%.2e", (a + b)}')
    done
    data6=("${total}")  # total wallclock
    util_write_file_with_header "${file_path}" "${append}"
    unset headers
    unset data0
}

parse_export_all_run_time_info_in_directory(){
    # parse run time info for all the cases in a directory
    # Inputs:
    # 	$1: directory
    #   $2: output file (optional)
    # File outout:
    # 	cases.log : log file of run time info
    local dir; local log_file;
    dir="$1"
    [[ -d "${dir}" ]] ||  cecho ${BAD} "${FUNCNAME[0]}: directory(${dir}) doesn't exist"
    [[ -n "$2" ]] && log_file="$2" || log_file="${dir}/cases.log"
    local is_case
    local if_second="0"
    for sub_dir in "${dir}"/*; do
    	if [[ -d "${sub_dir}" ]]; then
    		check_case "${sub_dir}" || parse_export_run_time_info "${sub_dir}" "${log_file}" "${if_second}"
    		[[ "${if_second}" = "0" ]] && if_second="1"  # in case this is the first entry
    	fi
    done
    echo "${FUNCNAME[0]}: log file generated: ${log_file}"
}

check_case(){
    # check a directory contains case data
    # Inputs:
    #	$1: directory
    # Return:
    # 	0: if this is an aspect case
    #   1: if not
    [[ -d "$1" ]] ||  cecho ${BAD} "${FUNCNAME[0]}: directory(${1}) doesn't exist"
    prm_path="$1/case.prm"
    [[ -e "${prm_file}" ]] && return 0 || return 1
}

check_job_id(){
    # check if a case is running on slurm
    # Inputs:
    #   $1: jobid
    local outputs=$(squeue -u "${USER}" -o %A)
    IFS=$'\n'; local job_ids=(${outputs})
    for job_id in ${job_ids[@]}; do
        if [[ ${job_id} =~ ^[0-9]+$ ]]; then 
            locate_workdir_with_id "$1"  # first locate where this is
            [[ $? == 1 ]] && return 1
            [[ -d "${return_value}" ]] ||  { cecho ${BAD} "${FUNCNAME[0]}: work directory(${return_value}) doesn't exist"; exit 1; }
            local path=$(fix_route "${return_value}")
            local case_path=$(fix_route "$1")
            printf "${path}\n${case_path}\n\n"
        fi
    done
}

check_case_running(){
    # check if a case is running
    # Inputs:
    #   $1: case path
    local case_path=$(fix_route "$1")
    local outputs=$(squeue -u "${USER}" -o %A)  # get all the job ids
    IFS=$'\n'; local job_ids=(${outputs})
    for job_id in ${job_ids[@]}; do
        if [[ ${job_id} =~ ^[0-9]+$ ]]; then 
            locate_workdir_with_id "${job_id}" 1  # first locate where this is
            [[ $? == 1 ]] && { cecho "$BAD" "${FUNCNAME[0]}: job id ${job_id} is not found in file system"; exit 1; }
            [[ -d "${return_value}" ]] ||  { cecho ${BAD} "${FUNCNAME[0]}: work directory(${return_value}) doesn't exist"; exit 1; }
            local path=$(fix_route "${return_value}")
	    [[ ${case_path} == ${path} ]] && { echo "case ${case_path} is running"; return 1; }  # case is running
	fi
    done
    return 0  # case not found
}

get_case_status(){
    # get the status of case (i.e. running (R), ended (E), terminate (T))
    # Inputs:
    #   $1: case_dir
    local case_dir="$1"
    prm_path="${case_dir}/case.prm"
    unset return_values
    util_get_prm_file_value "${prm_path}" "End time"
    local end_time="${return_values}"
    local status
    # [[  ]]
    echo "${status}"
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
    # $2 (0 or 1): quiet
    local outputs=$(scontrol show job "$1" | grep WorkDir)
    local quiet
    [[ -n "$2" ]] && quiet="$2" || quiet="0"
    if [[ -n ${outputs} ]]; then
        IFS='='; local entries=(${outputs})
	return_value=${entries[*]: -1}
	[[ ${quiet} == "0" ]] && printf "workdir: ${return_value}\n"
	return 0
    else
        cecho "${BAD}" "no job (id: $1) found"; return 1
    fi
}

################################################################################
# functions to submit cases
################################################################################

submit_case_peloton_p-billen(){
    ####
    # submit case to slurm
    # Inputs:
    #   $1: case directory
    now=$(pwd)
    cd "$1"
    eval "sbatch job_p-billen.sh"
    cd "${now}"
    return 0
}

submit_case_peloton_high2(){
    ####
    # submit case to slurm
    # Inputs:
    #   $1: case directory
    now=$(pwd)
    cd "$1"
    eval "sbatch job_high2.sh"
    cd "${now}"
    return 0
}

submit_case_stampede2_normal(){
    ####
    # submit case to slurm
    # Inputs:
    #   $1: case directory
    now=$(pwd)
    cd "$1"
    eval "sbatch job_normal.sh"
    cd "${now}"
    return 0
}

submit_case_stampede2_skx-normal(){
    ####
    # submit case to slurm
    # Inputs:
    #   $1: case directory
    now=$(pwd)
    cd "$1"
    eval "sbatch job_skx-normal.sh"
    cd "${now}"
    return 0
}

restart_case(){
    ####
    # restart a case
    # Inputs:
    #   $1: case directory
    #   $2: partition
    ####
    # change options in the case.prm file
    local case_dir="$1"
    local partition
    [[ -n "$2" ]] && partition="$2"
    local prm_file="${case_dir}/case.prm"
    [[ -e ${prm_file} ]] || { cecho "${BAD}" "No such file ${prm_file}"; exit 1; } 
    unset return_values  # rewrite prm file
    util_substitute_prm_file_contents "${prm_file}" "Resume computation" "true"
    if [[ ${partition} == "high2" ]]; then
        submit_case_peloton_high2 "${case_dir}"
    elif [[ ${partition} == "p-billen" ]]; then
        submit_case_peloton_p-billen "${case_dir}"
    else
        echo "no valid partition, skip restarting"
    fi
    # check for the restarted case
    return 0
}

# todo_restart
check_time_restart_case_combined(){
    return 0
}

check_time_restart_case(){
    ####
    # restart a case if the run time of that case is not reached
    # Inputs:
    #   $1: case directory
    #   $2: run time
    #   $3: partition
    ####
    local case_dir="$1"
    local time_plan="$2"
    local partition
    [[ -n "$3" ]] && partition="$3"
    # check if this is running
    check_case_running "$1" || return 0
    # read run_time
    local log_file="$1/output/log.txt"
    local restart_file="$1/output/restart.mesh"
    if [[ -e ${log_file} && -e ${restart_file} ]]; then 
        local outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_log_last_step" "${log_file}")
        local time=$(sed -E "s/^[^ ]*(\t|\ )*//g" <<< "${outputs}")
        if [[ $(eval "awk 'BEGIN{print ("${time}"<"${time_plan}")?1:0}'") -eq 1 ]]; then
	    outputs=$(awk -f "${ASPECT_LAB_DIR}/bash_scripts/awk_states/parse_snapshot" "${log_file}")
            IFS=$'\n'; local entries=(${outputs})
    	    IFS=' '; local return_value=(${entries[*]: -1})  # get the last line and split with ' '
	    [[ ${#return_value[*]} == 3 ]] || { cecho "${BAD}" "${FUNCNAME[0]}: entries from parsing snap shots mismatch, it's most likely there is old restart file in the folder, clean those up before proceed"; exit 1; }
            printf "Going to restart $1 at step ${return_value[0]} time ${return_value[1]}\n"
            restart_case "$1" "$3"
	else
	    printf "End time reached for $1 at ${time} (end time = ${time_plan})\n"
        fi
    else
        echo "${FUNCNAME[0]}: no snapshots to restart from, going to run from the start"
        if [[ ${partition} == "high2" ]]; then
            submit_case_peloton_high2 "${case_dir}"
    	elif [[ ${partition} == "p-billen" ]]; then
            submit_case_peloton_p-billen "${case_dir}"
    	elif [[ ${partition} == "normal" ]]; then
            submit_case_stampede2_normal "${case_dir}"
    	elif [[ ${partition} == "skx-normal" ]]; then
            submit_case_stampede2_skx-normal "${case_dir}"
        fi
    fi
    return 0
}

check_time_restart_case_in_directory(){
    ####
    # restart a case if the run time of that case is not reached
    # Inputs:
    #   $1: root directory
    #   $2: run time
    ####
    [[ -d "$1" ]] || { cecho "${BAD}" "no such directory $1"; exit 1; }
    [[ -n "$2" ]] || { cecho "${BAD}" "an end time must by assigned"; exit 1; }
    local partition
    [[ -n "$3" ]] && partition="$3"
    for sub_dir in "$1"/*; do
	local full_route=$(fix_route "${sub_dir}")
    	if [[ -d "${full_route}" ]]; then
    		check_case "${full_route}" ||  { echo "Find case ${full_route}"; check_time_restart_case "${full_route}" "$2" "${partition}"; }
        fi
    done
    return 0
}
    
copy_case_output_basics(){
    ####
    # copy case output, include only the basic entries and outputs
    # In this function, I copy the directory, keep the basename, and put it under a new target directory
    # Inputs:
    #   $1 - case_directory 
    #   $2 - target 
    ####
    [[ -d "$1" ]] || { cecho "${BAD}" "no such case directory $1"; exit 1; }
    local case_dir="$1"
    local base_name=$(basename "${case_dir}")
    [[ -d "$2" ]] || { cecho "${BAD}" "no such directory to output $2"; exit 1; }
    local target_dir="$2"

    local target_case_dir="$2/${base_name}"
    [[ -d "$target_case_dir" ]] || { mkdir "$target_case_dir"; }

    eval "cp ${case_dir}/*.prm ${target_case_dir}/" # copy the .prm file
    eval "cp ${case_dir}/*.wb ${target_case_dir}/" # copy the .wb file
    eval "cp ${case_dir}/*.sh ${target_case_dir}/" # copy the bash script files

    output_dir="${case_dir}/output"
    target_output_dir="${target_case_dir}/output"
    if [[ -d "${output_dir}" ]]; then
        [[ -d "$target_output_dir" ]] || { mkdir "$target_output_dir"; }
        eval "cp ${output_dir}/*.prm ${target_output_dir}/" # copy the .prm file
        eval "cp ${output_dir}/*.txt ${target_output_dir}/" # copy the .txt file
        eval "cp ${output_dir}/*.json ${target_output_dir}/" # copy the .json file
        eval "cp ${output_dir}/*.pvd ${target_output_dir}/" # copy the .pvd file
        eval "cp ${output_dir}/*.visit ${target_output_dir}/" # copy the .visit file
        eval "cp ${output_dir}/statistics ${target_output_dir}/" # copy the statistics file
    fi

    return 0
}

copy_case_output_by_vtu_snapshot(){
    ####
    # copy case output at a given vtu_snapshot
    # the vtu_snapshot here is the visualization step (step + initiation adaptive refinement level)
    # Inputs:
    #   $1 - case_directory 
    #   $2 - target 
    #   $3 - vtu_snapshot
    ####
    [[ -d "$1" ]] || { cecho "${BAD}" "no such case directory $1"; exit 1; }
    local case_dir="$1"
    [[ -d "$2" ]] || { cecho "${BAD}" "no such directory to output $3"; exit 1; }
    local target_dir="$2"
    local vtu_snapshot="$3"

    # first, only copy the basic components (e.g. .prm, .wb, statistics)
    copy_case_output_basics "${case_dir}" "${target_dir}"

    # then, make the solution directory and proceed by copying the files
    local output_dir="${case_dir}/output"
    local base_name=$(basename "${case_dir}")
    local solution_dir="${output_dir}/solution"
    local target_case_dir="${target_dir}/${base_name}"
    local target_output_dir="${target_case_dir}/output"
    local target_solution_dir="${target_output_dir}/solution"
    if [[ -d "${solution_dir}" ]]; then
        [[ -d "$target_solution_dir" ]] || { mkdir "$target_solution_dir"; }
        local stamp_vtu_snapshot=$(printf "%05d" ${vtu_snapshot})
        eval "cp ${solution_dir}/solution-${stamp_vtu_snapshot}.pvtu ${target_solution_dir}/"
        eval "cp ${solution_dir}/solution-${stamp_vtu_snapshot}.visit ${target_solution_dir}/"
        eval "cp ${solution_dir}/solution-${stamp_vtu_snapshot}.*.vtu ${target_solution_dir}/"
    fi
    
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
	[[ -n "$2" ]] || { cecho "${BAD}" "no case directory given (\$2)"; exit 1; }
	[[ -d "$2" ]] || { cecho "${BAD}" "directory given (\$2) doesn't exist"; exit 1; }
	log_file="$2/output/log.txt"
	[[ -e "${log_file}" ]] || { cecho "${BAD}" "log file (${log_file}) doesn't exist"; exit 1; }
	parse_run_time_info_last_step "${log_file}"
    elif [[ "$1" = "check_job_id" ]]; then 
        [[ -n "$2" ]] || { cecho "$BAD" "path of case (\$2) must be given"; exit 1; }
        check_job_id "$2"
    elif [[ "$1" = "check_case_running" ]]; then 
	check_case_running "$2"
    elif [[ "$1" = "case_info_with_id" ]]; then
    	[[ -n "$2" ]] || { cecho "${BAD}" "no id number given (\$2)"; exit 1; }
	    parse_run_time_info_last_step_with_id "$2"
    elif [[ "$1" = "export_case_info" ]]; then
    	[[ -n "$2" ]] || { cecho "${BAD}" "no case directory (\$2)"; exit 1; }
	    parse_export_run_time_info "$2" "$3"
    elif [[ "$1" = "export_all_case_info_in_directory" ]]; then
    	[[ -n "$2" ]] || { cecho "${BAD}" "no directory (\$2)"; exit 1; }
	parse_export_all_run_time_info_in_directory "$2"
    elif [[ "$1" = "all_case_info" ]]; then
	parse_all_time_info
    elif [[ "$1" = "case_info_in_dir" ]]; then
	parse_run_time_info_last_step_in_dir "$2"
    elif [[ "$1" == "restart" ]]; then
        restart_case "$2"
    elif [[ "$1" == "check_time_restart" ]]; then
        [[ -n "$2" && -n "$3" ]] || \
        { cecho "$BAD" "path of case (\$2) and time (\$3) must be given"; exit 1; }
        check_time_restart_case "$2" "$3" "$4"
    elif [[ "$1" == "check_time_restart_directory" ]]; then
        check_time_restart_case_in_directory "$2" "$3" "$4"
    elif [[ "$1" == "copy_case_output_by_vtu_snapshot" ]]; then
        copy_case_output_by_vtu_snapshot "$2" "$3" "$4"
    else
    	cecho "${BAD}" "option ${1} is not valid\n"
    fi
}

set +a  # return to default setting

if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
