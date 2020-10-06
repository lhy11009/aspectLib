#!/bin/bash
# case manager
# Usage:
#   ./aspect_lib.sh + command + options
# future: use a file to compile remote address

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"

source "${dir}/utilities.sh"
source "${dir}/extra_configurations.sh"


################################################################################
# parse parameters from command line
# Inputs:
#   $1: options
#        formate of options should be "-a val1 --b=valb ..."
parse_options(){

    # parse options
    while [ -n "$1" ]; do
      param="$1"
      case $param in
        -h|--help)
          usage  # help information
          exit 0
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
        # list
        #####################################
        -l)
          shift
          vlist=()
          while [[ true ]]; do
            [[ -z "$1" ]] && { cecho ${BAD} "${FUNCNAME[0]}: no value given for list"; exit 1; }
            vlist+=("$1")
            [[ -z "$2" || "$2" = -* ]] && break || shift
          done
        ;;
      esac
      shift
    done

    # check values
    [[ -z ${bool} || ${bool} = "true" || ${bool} = "false" ]] || { cecho ${BAD} "${FUNCNAME[0]}: bool value must be true or false"; exit 1; }
}


################################################################################
# help message
usage()
{
    printf "\
Submit a job to cluster with a slurm system

Usage:
  %s [options] [server_name] [file_name]

Options:
    --bool -b       A bool value that is either 'true' or 'false'
    
    -l              a lisl of value
                    for example:
                        -l 1.0 2.0 3.0
"
}



create_case(){
    # create a case locally
    local py_script="$1"
    local local_root="$2"
    eval "python -m ${py_script} create -j config_case.json >.temp 2>&1"  # debug
    [[ $? -eq 1 ]] && { cecho ${BAD} "${py_script} failed"; exit 1; }
    # get case name
    _info=$(cat ".temp")
    case_name=$(echo "${_info}" | sed -n '2'p)
    # assertion
    local case_dir="${local_root}/${case_name}"
    local case_prm="${case_dir}/case.prm"
    [[ -d ${case_dir} && -e ${case_prm} ]] || { cecho ${BAD} "Case generation failed"; exit 1; }
    cecho ${GOOD} "${_info}"
}


create_group(){
    # create a group of case locally
    local py_script="$1"
    local local_root="$2"
    eval "python -m ${py_script} create_group -j config_group.json 2>&1 > .temp"
    [[ $? -eq 1 ]] && { cecho ${BAD} "${py_script} failed"; exit 1; }
    # get case names
    _info=$(cat ".temp")
    local group_name=$(echo "${_info}" | sed -n '2'p)
    IFS=' ' read -r -a case_names <<< $(cat ".temp" | sed -n '1,2!'p)
    group_dir="${local_root}/${group_name}"
    for case_name in ${case_names[@]}; do
        local case_dir="${group_dir}/${case_name}"
        local case_prm="${case_dir}/case.prm"
        [[ -d ${case_dir} && -e ${case_prm} ]] || { cecho ${BAD} "Case generation failed"; exit 1; }
    done
    cecho ${GOOD} "${_info}"
}


submit(){
    # future parse from json file
    local case_dir="$1"
    local case_name=$(basename "${case_dir}")
    local romote_case_dir="$2"
    local server_info="$3"
    local flag=''  # a vacant flag for adding optional parameters
    local case_prm="${case_dir}/case.prm"
    local remote_case_prm="${remote_case_dir}/case.prm"
    # output machine time output
    local remote_time_file="${remote_case_dir}/output/machine_time"
    flag="${flag} -lt ${remote_time_file}"

    # get configuration from a file
    local node=1  # backward compatible
    total_tasks=$(sed -n '1'p "slurm_config")
    time_by_hour=$(sed -n '2'p "slurm_config")
    partition=$(sed -n '3'p "slurm_config")
    nodes=$(sed -n '4'p "slurm_config")

    # scp to remote
    local remote_target=$(dirname "${remote_case_dir}")
    eval "${RSYNC} -r ${case_dir} ${server_info}:${remote_target}"

    # check file arrival
    local status_
    while [[ true ]]; do
        ssh ${server_info} << EOF > '.temp'
            eval "[[ -e ${remote_case_prm} ]] && echo \"0\" || echo \"1\""
EOF
        status_=$(cat '.temp'| sed -n '$'p)
        ((status_==0)) && break || { cecho ${WARN} "Files haven't arrived yet, sleep for 2s"; sleep 2s; }
    done

    # add an optional log file
    [[ "$4" != '' ]] && flag="${flag} -l $4"  # add -l log_file to flag, if $4 given
    # submit using slurm.sh,
    # determine if there is a valid job id, future
    # also add -P option for project name
    ssh ${server_info} << EOF > '.temp'
        eval "slurm.sh -N ${nodes} -n ${total_tasks} -t ${time_by_hour} -p ${partition} -P ${project} ${remote_case_prm} ${flag}"
EOF
    # get job_id
    local _info=$(cat '.temp'| sed -n '$'p)
    local job_id=$(echo "${_info}" | sed 's/Submitted\ batch\ job\ //')
    if ! [[ ${job_id} != '' && ${job_id} =~ ^[0-9]*$  ]]; then
        cecho ${BAD} "submit case: ${case_name} failed"
        return 1
    else
        cecho ${GOOD} "submit case: ${case_name} succeeded, job id: ${job_id}"
        echo "${job_id}" > ".temp"  # use .temp to transfer information
        return 0
    fi
}


install(){
    local project="$1"
    local server_info="$2"

    # new folder
    local project_dir="${ASPECT_PROJECT_DIR}/${project}"
    [[ -d "${project_dir}" ]]  || mkdir "${project_dir}"

    # on server side
    get_remote_environment "${server_info}" "ASPECT_PROJECT_DIR"
    local remote_project_dir="${return_value}/${project}"
    ssh ${server_info} << EOF
        eval "[[ -d \"${remote_project_dir}\" ]]  || mkdir \"${remote_project_dir}\" "
EOF

    # set alias, add a line every time it executes, future: fix this bug
    echo "export ${project}_DIR=\"${project_dir}\"" >> "${dir}/env/enable.sh"
    echo "export ${project}_DIR=\"${remote_project_dir}\"" >> "${dir}/env/enable_peloton.sh";

    # new mkdocs project
    previous=$(pwd)
    cd "${project_dir}"
    if ! [[ -d "mkdocs_project" ]]; then
        eval "mkdocs new mkdocs_project"
        yml_file="${project_dir}/mkdocs_project/mkdocs.yml"  # make a .yml file
        echo "site_name: ${project}" > "${yml_file}"
        echo "nav:" >> "${yml_file}"
        echo "    - Home: index.md" >> "${yml_file}"
        echo "theme: readthedocs" >> "${yml_file}"
    fi
    cd "${previous}"
}

################################################################################
# Translate a visit script
# Inputs:
#    filein: name of the script
#    keys: keys to translate
#    values: values to translate
tranlate_visit_script(){
    # check variable
    # check_variable 'keys'
    # check_variable 'values'
    # check_variable 'filein'
    # check_variable 'fileout'

    # read file
    contents=$(cat ${filein})

    # do substutions
    local i=0
    for key in ${keys[@]}; do
        value=${values[i]}
        contents=${contents//$key/$value}  # substitute key with value
        ((i++))
    done

    # output
    echo "${contents}" > "${fileout}"
}


################################################################################
# Translate a visit script, run it and generate plots
# Inputs:
#    fileins: list name of the script
#    file_keys_values: file that holds keys and values to substitute
plot_visit_scripts(){
    # future
    echo '0'
}
        

################################################################################
# Run tests
# Inputs:
run_tests(){
    current_dir=$(pwd)

    # run python tests
    cd ${ASPECT_LAB_DIR}
    # eval "python -m pytest tests -v --cov"

    cd ${current_dir}

    # run bash tests
    bash_tests_dir="${ASPECT_LAB_DIR}/bash_tests"
    for file in ${bash_tests_dir}/* ; do
        [[ ${file} =~ /test.*\.sh$ && -x ${file} ]] && eval "${file} ${project} ${server_info}"
    done
}


################################################################################
# Generate visit plots for a single case
# Inputs
#    case_dir: case directory
#    bool: whether to plot
#       if this option is 'false', this function will just translate the scripts without plotting
plot_visit_case(){
    # check folders
    [ -d "${case_dir}" ] || { cecho ${BAD} "plot_visit_case: Case folder - ${case_dir} doesn't exist"; exit 1; }
    data_sub_dir="${case_dir}/output"
    [ -d "${data_sub_dir}" ] || { cecho ${BAD} "plot_visit_case: Data folder - ${data_sub_dir} doesn't exist"; exit 1; }
    # dir for transfered visit scripts
    local visit_temp_dir="${dir}/visit_scripts_temp"
    [ -d "${visit_temp_dir}" ] || mkdir "${visit_temp_dir}" ]
    # dir for image output
    local img_dir="${case_dir}/img"
    [ -d "${img_dir}" ] || mkdir "${img_dir}" ]

    # get a list of scripts to plot
    visit_script_bases=("initial_slab.py" "export_particles.py" "slab.py")
    visit_script_dir="${dir}/visit_scripts/${project}"
        
    # call python module to generate visit_keys_values file
    eval "python -m shilofue.${project} bash_options -i ${case_dir}"
    
    # get keys and values
    keys_values_file="${dir}/visit_keys_values"
    [ -r "${keys_values_file}" ] || { cecho ${BAD} "plot_visit_case: Files containing keys and values - ${keys_values_file} cannot be read"; exit 1; }
    read_keys_values "visit_keys_values"
    
    # do substitution and run
    for visit_script_base in ${visit_script_bases[@]}; do
        filein="${visit_script_dir}/${visit_script_base}"
        fileout="${visit_temp_dir}/${visit_script_base}"

        # translate script
        tranlate_visit_script

        # run
        [[ -z ${bool} || ${bool} = "true" ]] && echo "exit()" | eval "visit -nowin -cli -s ${fileout}"
    done
}


################################################################################
# outputs from solver
# Inputs:
#   filein: file to read from
#   fileout: file to output to
#   vlist: [start step, interval, end step].
#       If this is none, take all steps   
parse_solver_output(){
    # parse in options
    # assert the list of value
    local start; local interval; local end;
    if [[ -n ${vlist} ]]; then
        if [[ ${#vlist[@]} -eq 3 ]]; then
            start=${vlist[0]}
            interval=${vlist[1]}
            end=${vlist[2]}
        else
            cecho ${BAD} "${FUNCNAME[0]}: If vlist exists, it must be a list of length of 3"; exit 1;
        fi
    else
        start=0
        interval=1
        # set a max number to avoid dead loop
        end=1000000
    fi

    # loop for timesteps
    timestep=${start}
    while ((timestep<end)); do
        parse_solver_output_timestep
        
        # check if this is successful
        # we take for granted that any non-zero state
        # indicates that we reach the EOF
        [[ $? -eq 0 ]] || break

        ((timestep+=interval))
    done
    return 0
}


################################################################################
# outputs from solver by giving a case directory
# Before converting result, we check if there is already results presented and we
# check the relative time of that file to the new output file.
# Inputs:
#   case_dir: directory to parse from
#   vlist: [start step, interval, end step].
#       If this is none, take all steps   
# Returns:
#   0: normal
#   -1: no operation needed
parse_case_solver_output()
{
    # unset
    unset filein
    
    # find newest .stdout file
    local temp; local id; local filename
    local idin=0
    for file_ in "${case_dir}"/*"stdout"; do
        # fix bug
        [[ -e "${file_}" ]] || break
        # get id
        filename=$(basename "${file_}")
        temp=${filename#*"-"}
        id=${temp%%".stdout"}
        ((id>idin)) && { filein="${file_}"; ((idin=id)); }
    done
    [[ -z "${filein}" ]] && { cecho ${WARN} "${FUNCNAME[0]}: fail to get file in dir ${case_dir}"; return 1; }

    # check output dir 
    local output_dir="${case_dir}/output"
    [[ -d ${output_dir} ]] || mkdir "${output_dir}"

    # check results existence and compare time
    fileout="${output_dir}/solver_output"
    [[ -e "${fileout}" && "${fileout}" -nt ${filein} ]] && return -1

    # call parse_solver_output function
    parse_solver_output
    return 0
}


################################################################################
# outputs from solver for a timestep
# Inputs:
#   filein: file to read from
#   fileout: file to output to
#   timestep(int): time step
# Returns:
#   0
#   1: reach end of the file
#   2: no such outputs presented
parse_solver_output_timestep(){
    # read input 
    [[ -z ${filein} ]] && { cecho ${BAD} "${FUNCNAME[0]}: filein must be given"; exit 1; }
    [[ -z ${fileout} ]] && { cecho ${BAD} "${FUNCNAME[0]}: fileout must be given"; exit 1; }
    [[ -z ${timestep} ]] && { cecho ${BAD} "${FUNCNAME[0]}: timestep must be given"; exit 1; }
    
    # parse content from stdout file 
    parse_stdout1 "${filein}" "${timestep}"
    [[ $? -eq 0 ]] || { echo "${FUNCNAME[0]}: timestep ${timestep} seems to hit end of the file ${filein}"; return 1; }

    # parse one time step
    local solver_outputs=()
    local output=''
    local start=0
    parse_block_outputs "${content}" "Rebuilding Stokes"

    # get useful information
    # number of nonlinear solver
    local nnl="${#block_outputs[@]}"
    # Relative nonlinear residual (total Newton system)
    local rnrs=()
    # norm of the rhs
    local nors=()
    # newton_derivative_scaling_factor
    local ndsfs=()
    local temp; local temp1
    for block_output in "${block_outputs[@]}"; do
        # get rnr
        parse_output_value "${block_output}" "nonlinear iteration" ":" ","
        [[ -z ${value} ]] && { cecho ${WARN} "${FUNCNAME[0]}: ${filein} doens't have solver outputs"; return 2; } 
        rnrs+=("${value}")
        # get nors
        parse_output_value "${block_output}" "norm of the rhs" ":" ","
        nors+=("${value}")
        # get ndrsf
        parse_output_value "${block_output}" "newton_derivative_scaling_factor" ":" ","
        # if value is not present, append by 0
        [[ -n ${value} ]] && ndsfs+=("${value}") || ndsfs+=("0")
    done

    # output header if file doesn't exist
    if ! [[ -e ${fileout} ]]; then 
        printf "# 1: Time step number\n" >> "${fileout}"
        printf "# 2: Index of nonlinear iteration\n" >> "${fileout}"
        printf "# 3: Relative nonlinear residual\n" >> "${fileout}"
        printf "# 4: Norms of the rhs\n" >> "${fileout}"
        printf "# 5: Newton Derivative Scaling Factor\n" >> "${fileout}"
    fi

    # get length of array
    local length=${#rnrs[@]}
 
    # output
    local i=0
    while ((i<length)); do
        # output to file 
        printf "%-15s %-15s %-15s %-15s %s\n" "${timestep}" "${i}" "${rnrs[$i]}" "${nors[$i]}" "${ndsfs[$i]}" >> "${fileout}"
        ((i++))
    done
    
    return 0
}


################################################################################
# future
# build a project in aspect
# usage of this is to bind up source code and plugins
# Inputs:
#   project: name of the project
#   $1: release or debug, default is debug following aspect's routine
build_aspect_project(){
    build_dir="${ASPECT_SOURCE_DIR}/build_${project}"
    [[ -d ${build_dir} ]] || mkdir ${build_dir}
    local mode
    if [[ -n $1 ]]; then
        [[ $1="debug" || $1="release" ]] || { cecho ${BAD} "${FUNCNAME[0]}: mode is either \'debug\' or \'release\'"; exit 1; }
        mode=$1
    else
        mode="debug"
    fi

    # get the project json file
    json="${ASPECT_LAB_DIR}/files/${project}/project.json"
    [[ -e ${json} ]] || cecho ${WARN} "${FUNCNAME[0]}: json file of project(i.e. ${json}) doesn't exist"

    # get the list of plugins
    plugins=("prescribe_field" "subduction_temperature2d")

    # copy plugins
    plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    for plugin in ${plugins[@]}; do
        # check plugin existence
        plugin_dir="${plugins_dir}/${plugin}"
        [[ -d ${plugin_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: plugin(i.e. ${plugin_dir}) doesn't exist"; exit 1; }

        # remove old ones
        plugin_to_dir="${build_dir}/${plugin}"
        [[ -d ${plugin_to_dir} ]] && rm -r ${plugin_to_dir}

        # copy new ones
        eval "cp -r ${plugin_dir} ${build_dir}/"
        cecho ${GOOD} "${FUNCNAME[0]}: copyied plugin(i.e. ${plugin})"
    done

    # build
    local current_dir=$(pwd)
    # Here we pick nproc - 1, this make sure that we don't use up all resources. 
    # But this will cause problem when nproc = 1
    local nproc=$(($(nproc)-1))
    cd ${build_dir}
    # build source code
    eval "cmake .."
    quit_if_fail "${FUNCNAME[0]}: cmake inside ${build_dir} failed"
    eval "make ${mode}"
    quit_if_fail "${FUNCNAME[0]}: \"make ${mode}\" inside ${build_dir} failed"
    eval "make -j ${nproc}"
    quit_if_fail "${FUNCNAME[0]}: make inside ${build_dir} failed"
    # build plugins
    for plugin in ${plugins[@]}; do
        plugin_to_dir="${build_dir}/${plugin}"
        cd ${plugin_to_dir}
        # remove cache before compling
        [[ -e "${plugin_to_dir}/CMakeCache.txt" ]] && eval "rm ${plugin_to_dir}/CMakeCache.txt"
        eval "cmake -DAspect_DIR=${build_dir}"
        quit_if_fail "${FUNCNAME[0]}: cmake inside ${plugin_to_dir} failed"
        eval "make"
        quit_if_fail "${FUNCNAME[0]}: make inside ${plugin_to_dir} failed"
    done
}


main(){
    # parameter list, future
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local py_script="shilofue.${project}"

    # parse commend
    _commend="$2"

    # parse options
    parse_options $@

    # check project
    [[ -d ${local_root} || ${_commend} = 'install' || ${_commend} = 'write_time_log' || ${_commend} = 'keep_write_time_log' ]] || { cecho ${BAD} "Project ${project} is not included"; exit 1; }

    # execute
    if [[ ${_commend} = 'install' ]]; then
        [[ "$#" -eq 3 ]] || { cecho ${BAD} "for install, server_info must be given"; exit 1; }
        local server_info="$3"
        install "${project}" ${server_info}

    elif [[ ${_commend} = 'create' ]]; then
        create_case "${py_script}" "${local_root}"

    elif [[ ${_commend} = 'create_group' ]]; then
        create_group "${py_script}" "${local_root}"

    elif [[ ${_commend} = 'submit' ]]; then
        local case_name="$3"
        local server_info="$4"
        local case_dir="${local_root}/${case_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"} # substitution
        local log_file="$5"  # add an optional log_file
        if [[ "${log_file}" != '' ]]; then
            # if there is no $5 given, log file is ''
            log_file=$(fix_route "${log_file}")
            log_file=${log_file/"${local_root}"/"${remote_root}"} # substitution
        fi
        submit "${case_dir}" "${remote_case_dir}" "${server_info}" "${log_file}"
        quit_if_fail "aspect_lib.sh submit failed for case ${case_name}"

    elif [[ ${_commend} = 'submit_group' ]]; then
        local group_name="$3"
        local server_info="$4"
        local group_dir="${local_root}/${group_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        local remote_group_dir=${group_dir/"${local_root}"/"${remote_root}"}
        ssh "${server_info}" eval "[[ -d ${remote_group_dir} ]] && { rm -r ${remote_group_dir}; mkdir ${remote_group_dir}; }|| mkdir ${remote_group_dir}"
        local log_file="$5"  # add an optional log_file, future, move this to global settings
        if [[ "${log_file}" != '' ]]; then
            # if there is no $5 given, log file is ''
            log_file=$(fix_route "${log_file}")
            log_file=${log_file/"${local_root}"/"${remote_root}"} # substitution
        fi
        # get a list of cases and submit
        local job_ids=""
        for case_dir in "${group_dir}/"*; do
            if [[ -d "${case_dir}" ]]; then
                # select directories
                local _files=$(ls "${case_dir}")
                if [[ "${_files[@]}" =~ 'case.prm' ]]; then
                    local remote_case_dir=${case_dir/"${local_root}"/"${remote_root}"}
                    # call submit functions
                    submit "${case_dir}" "${remote_case_dir}" "${server_info}" "${log_file}"
                    quit_if_fail "aspect_lib.sh submit group failed for case ${case_dir}"
                    local job_id=$(cat ".temp")  # get job id
                    job_ids="${job_ids} ${job_id}"
                fi
            fi
        done
        echo "${job_ids[@]}" > ".temp"
        return 0

    elif [[ ${_commend} = 'create_submit' ]]; then
        local server_info="$3"
        local log_file="$4"  # optional log file
        ./aspect_lib.sh "${project}" 'create'
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh create failed"; exit 1; }
        # get case name
        local _info=$(cat ".temp")
        local case_name=$(echo "${_info}" | sed -n '2'p)
        # submit to server
        ./aspect_lib.sh "${project}" 'submit' "${case_name}" "${server_info}" "${log_file}"
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh submit failed"; exit 1; }

    elif [[ ${_commend} = 'create_submit_group' ]]; then
        local server_info="$3"
        local log_file="$4"  # optional log file
        ./aspect_lib.sh "${project}" 'create_group'
        [[ $? -eq 0 ]] || {  cecho ${BAD} "aspect_lib.sh create failed"; exit 1; }
        # get group name
        local _info=$(cat ".temp")
        local group_name=$(echo "${_info}" | sed -n '2'p)
        # call self
        ./aspect_lib.sh "${project}" 'submit_group' "${group_name}" "${server_info}" "${log_file}"
        quit_if_fail "aspect_lib.sh submit_group failed"

    elif [[ ${_commend} = 'terminate' ]]; then
        # future
        echo '0'

    elif [[ ${_commend} = 'remove' ]]; then
        # future
        echo '0'

    elif [[ ${_commend} = 'translate_visit' ]]; then
        # translate visit scripts
        filein="$3"
        filein_base=$(basename ${filein})
        fileout="${ASPECT_LAB_DIR}/visit_scripts_temp/${filein_base}"

        # get keys
        read_keys_values "${ASPECT_LAB_DIR}/visit_keys_values"

        # call function
        tranlate_visit_script
    
    elif [[ ${_commend} = 'plot_visit_case' ]]; then
        # plot visit for a case
        # example command line:
        # ./aspect_lib.sh TwoDSubduction plot_visit_case $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12
        # no plot so that we could run in gui later:
        #   ./aspect_lib.sh TwoDSubduction plot_visit_case $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12 -b false
        case_dir="$3"
        [[ -d ${case_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: case directory(i.e. ${case_dir}) doesn't exist"; exit 1; }

        # call function
        plot_visit_case

    elif [[ ${_commend} = 'parse_solver_output' ]]; then
        # parse solver information from stdout file
        # example command line:
        # ./aspect_lib.sh TwoDSubduction parse_solver_output $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12/task.stdout ./solver_output
        # only output first 20 steps:
        # ./aspect_lib.sh TwoDSubduction parse_solver_output $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12/task.stdout ./solver_output -l 0 1 20
        filein="$3"
        fileout="$4"
        
        # call function
        parse_solver_output
    
    elif [[ ${_commend} = 'parse_case_solver_output' ]]; then
        # parse solver information from stdout file by giving a case directory
        # the .stdout file must be placed under this directory
        # and the output file 'solver_output' goes into the 'output' directory
        # example command line:
        # ./aspect_lib.sh TwoDSubduction parse_case_solver_output $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12
        # only output first 20 steps:
        # ./aspect_lib.sh TwoDSubduction parse_case_solver_output $TwoDSubduction_DIR/isosurf_global2/isosurfULV3.000e+01testS12 -l 0 1 20
        case_dir="$3"
        
        # call function
        parse_case_solver_output

    elif [[ ${_commend} = 'write_time_log' ]]; then
        # write time and machine time output to a file
        # example command line:
        # ./aspect_lib.sh foo write_time_log /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS13\
        # 2537585 /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS13/output/machine_time
        [[ -n $3 && -d $3 ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$3 must be a valid directory"; exit 1; }
        [[ -n $4 && $4=~^[0-9]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$4 must be a valid job id"; exit 1; }
        [[ -n $5 ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$5 must be a valid path of a file"; exit 1; }
        write_time_log $3 $4 $5
    
    
    elif [[ ${_commend} = 'keep_write_time_log' ]]; then
        # write time and machine time output to a file
        # example command line:
        # nohup ./aspect_lib.sh foo keep_write_time_log /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS13\
        # 2537585 /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS13/output/machine_time &
        [[ -n $3 && -d $3 ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$3 must be a valid directory"; exit 1; }
        [[ -n $4 && $4=~^[0-9]+$ ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$4 must be a valid job id"; exit 1; }
        [[ -n $5 ]] || { cecho ${BAD} "${FUNCNAME[0]}: write_time_log, \$5 must be a valid path of a file"; exit 1; }

        while true
        do
            write_time_log $3 $4 $5
            [[ $? -eq 0 ]] || { printf "${FUNCNAME[0]}: stop writing time log\n"; exit 0; }
            sleep 1h
        done
    
    elif [[ ${_commend} = 'build' ]]; then
        # build a project in aspect
        # usage of this is to bind up source code and plugins
        # example command line:
        #   local build:
        #       ./aspect_lib.sh TwoDSubduction build
        #       ./aspect_lib.sh TwoDSubduction build debug
        #       ./aspect_lib.sh TwoDSubduction build release
        build_aspect_project $3
    
    elif [[ ${_commend} = 'build_remote' ]]; then
        #   server build(add a server_info):
        #   example command lines:
        #       ./aspect_lib.sh TwoDSubduction build lochy@peloton.cse.ucdavis.edu
        #       ./aspect_lib.sh TwoDSubduction build lochy@peloton.cse.ucdavis.edu debug
        #       ./aspect_lib.sh TwoDSubduction build lochy@peloton.cse.ucdavis.edu release
        server_info="$3"
        ssh ${server_info} << EOF
            eval "\${ASPECT_LAB_DIR}/aspect_lib.sh ${project} build $4"
EOF

    elif [[ ${_commend} = 'test' ]]; then
        # run tests
        # example command line:
        #   local test:
        #       ./aspect_lib.sh TwoDSubduction test
        #   server test:
        #       ./aspect_lib.sh TwoDSubduction test lochy@peloton.cse.ucdavis.edu
        # get server info
        server_info="$3"

        # call scripts in bash_tests folder
        run_tests
    
    else
        cecho ${BAD} "Bad commend: ${_commend}"
    fi
    return 0
}

main $@
