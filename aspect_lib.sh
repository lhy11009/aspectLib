#!/bin/bash
# case manager
# Usage:
#   ./aspect_lib.sh + command + options
# todo_future: use a file to compile remote address

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"

source "${dir}/utilities.sh"
source "${dir}/extra_configurations.sh"

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
    local case_dir="$1"
    local case_name=$(basename "${case_dir}")
    local romote_case_dir="$2"
    local server_info="$3"
    local flag=''  # a vacant flag for adding optional parameters
    local case_prm="${case_dir}/case.prm"
    local remote_case_prm="${remote_case_dir}/case.prm"
    # get configuration from a file
    total_tasks=$(sed -n '1'p "slurm_config")
    time_by_hour=$(sed -n '2'p "slurm_config")
    partition=$(sed -n '3'p "slurm_config")
    # scp to remote
    local remote_target=$(dirname "${remote_case_dir}")
    eval "${RSYNC} -r ${case_dir} ${server_info}:${remote_target}"
    sleep 25s
    # add an optional log file
    [[ "$4" != '' ]] && flag="${flag} -l $4"  # add -l log_file to flag, if $4 given
    # submit using slurm.sh,
    # determine if there is a valid job id, todo
    # also add -P option for project name
    ssh ${server_info} << EOF > '.temp'
        eval "slurm.sh -n ${total_tasks} -t ${time_by_hour} -p ${partition} -P ${project} ${remote_case_prm} ${flag}"
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

    # set alias, add a line every time it executes, todo_future: fix this bug
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
    printf "${contents}" > "${fileout}"
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
    visit_script_bases=("initial_slab.py" "export_particles.py")
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
        echo "exit()" | eval "visit -nowin -cli -s ${fileout}"
    done
}


################################################################################
# todo
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
    # parameter list, todo
    local project="$1"
    local local_root=$(eval "echo \${${project}_DIR}")
    local py_script="shilofue.${project}"
    _commend="$2"
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
        # todo
        local group_name="$3"
        local server_info="$4"
        local group_dir="${local_root}/${group_name}"
        # get remote case directory
        get_remote_environment "${server_info}" "${project}_DIR"
        local remote_root=${return_value}
        local remote_group_dir=${group_dir/"${local_root}"/"${remote_root}"}
        ssh "${server_info}" eval "[[ -d ${remote_group_dir} ]] && { rm -r ${remote_group_dir}; mkdir ${remote_group_dir}; }|| mkdir ${remote_group_dir}"
        local log_file="$5"  # add an optional log_file, todo_future, move this to global settings
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
        case_dir="/home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12"

        # call function
        plot_visit_case


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
