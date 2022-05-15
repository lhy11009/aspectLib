#!/bin/bash


################################################################################
# Build aspect branch
#
# Dependencies:
#    utilities.sh
#    env:
#        ASPECT_SOURCE_DIR
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"


usage(){
  # usage of this script
    _text="
${BASH_SOURCE[0]}

(Descriptions)

Dependencies:
   env:
       ASPECT_LAB_DIR
       ASPECT_SOURCE_DIR

Example Usage:

    build all (all plugins, release and debug mode):

        Lib_build_aspect all master_TwoD

        the last one is the branch to build

     build a plugin:
        (build it in a build directory of aspect)
        Lib_build_aspect subduction_temperature2d TwoDSubduction
    
        (build it in a seperate folder)
        Lib_build_aspect visco_plastic_TwoD master_TwoD ${ASPECT_SOURCE_DIR}/build_plugins

    build all plugins separately
        Lib_build_aspect all_plugins master_TwoD_hefesto

    fix test
        Lib_build_aspect fix_test "${ASPECT_SOURCE_DIR}" "build_master_TwoD" TwoDSubduction_HeFESTo_steinberg
    


"
    printf "${_text}"

}


build_aspect_project(){
    ###
    # build a project in aspect
    # usage of this is to bind up source code and plugins
    # Inputs:
    #   project: name of the project
    ###
    local build_dir
    [[ $1 == "master" ]] && build_dir="${ASPECT_SOURCE_DIR}/build" || build_dir="${ASPECT_SOURCE_DIR}/build_$1"
    local mode
    if [[ -n $2 ]]; then
        [[ $2 == "debug" || $2 == "release" ]] || { cecho ${BAD} "${FUNCNAME[0]}: mode is either \'debug\' or \'release\'"; exit 1; }
        mode="$2"
    else
        mode="debug"
    fi
    [[ ${mode} == "debug" ]] && build_dir="${build_dir}_debug"
    [[ -d ${build_dir} ]] || mkdir ${build_dir}

    # shift the git dir
    cd "${ASPECT_SOURCE_DIR}"
    echo "checkout $1 branch"
    eval "git checkout $1"
    quit_if_fail "${FUNCNAME[0]}: git checkout doesn't work"
    cd "${ASPECT_LAB_DIR}"

    # get the list of plugins
    local plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    # [[ -d ${plugins_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ${plugins_dir} doesn't exist"; exit 1; }
    plugins=()
    for folder in ${plugins_dir}/*; do
        if [[ -d ${folder} ]]; then
            plugin=$(basename "$folder")
            plugins+=("${plugin}")
        fi
    done
    echo "plugins: ${plugins[@]}"  # debug

    # build
    local current_dir=$(pwd)
    # Here we pick nproc - 1, this make sure that we don't use up all resources.
    # But this will cause problem when nproc = 1
    # local nproc=$(nproc)
    local nproc=$(nproc)
    (( ${nproc} > 6 )) && nproc=6
    # local nproc=8
    cd ${build_dir}
    # build source code
    local cmake_appendix=""
    [[ -n $WORLD_BUILDER_SOURCE_DIR ]] && cmake_appendix="${cmake_appendix} -DWORLD_BUILDER_SOURCE_DIR=$WORLD_BUILDER_SOURCE_DIR"
    echo "cmake .. ${cmake_appendix}"
    eval "cmake .. ${cmake_appendix}"
    quit_if_fail "${FUNCNAME[0]}: cmake inside ${build_dir} failed"
    echo "make ${mode}"
    eval "make ${mode}"
    quit_if_fail "${FUNCNAME[0]}: \"make ${mode}\" inside ${build_dir} failed"
    echo "make -j ${nproc}"  # screen output
    eval "make -j ${nproc}"
    quit_if_fail "${FUNCNAME[0]}: make inside ${build_dir} failed"

    # build plugins
    for plugin in ${plugins[@]}; do
        build_aspect_plugin "${build_dir}" "${plugin}"
    done
}

build_all_plugins_separately(){
    # todo
    # build plugins separately for a project in aspect
    # usage of this is to bind up source code and plugins
    # built objects are aimed to be used in tests
    # Inputs:
    #   project: name of the project
    local build_dir="${ASPECT_SOURCE_DIR}/build_$1"
    local build_plugin_dir="${ASPECT_SOURCE_DIR}/build_plugins"
    local plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    [[ -d ${plugins_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ${plugins_dir} doesn't exist"; exit 1; }
    [[ -d ${build_plugin_dir} ]] || mkdir "${build_plugin_dir}"
    # list of plugins
    plugins=()
    for folder in ${plugins_dir}/*; do
        if [[ -d ${folder} ]]; then
            plugin=$(basename "$folder")
            plugins+=("${plugin}")
        fi
    done
    echo "plugins: ${plugins[@]}"  # debug
    # build
    for plugin in ${plugins[@]}; do
        build_aspect_plugin "${build_dir}" "${plugin}" "${build_plugin_dir}"
    done
}

build_aspect_plugin(){
    ###
    # build a plugin with source code in aspect
    # usage of this is to bind up source code and plugins
    # Inputs:
    #   project: name of the project
    #   $1: name of plugin
    ###
    local build_dir="$1"
    [[ -d ${build_dir} ]] || mkdir ${build_dir}
    local plugin="$2"
    local target_dir
    [[ -z $3 ]] && target_dir="${build_dir}" || target_dir="$3"

    # copy plugins
    plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    # check plugin existence
    plugin_dir="${plugins_dir}/${plugin}"
    [[ -d ${plugin_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: plugin(i.e. ${plugin_dir}) doesn't exist"; exit 1; }

    # remove old ones
    plugin_to_dir="${target_dir}/${plugin}"
    echo "plugin_to_dir: ${plugin_to_dir}"  # debug
    [[ -d ${plugin_to_dir} ]] && rm -r ${plugin_to_dir}
    # copy new ones
    eval "cp -r ${plugin_dir} ${target_dir}/"
    cecho ${GOOD} "${FUNCNAME[0]}: copyied plugin(i.e. ${plugin})"

    # build
    cd ${plugin_to_dir}
    # remove cache before compling
    [[ -e "${plugin_to_dir}/CMakeCache.txt" ]] && eval "rm ${plugin_to_dir}/CMakeCache.txt"
    echo "cmake -DAspect_DIR=${build_dir}"
    eval "cmake -DAspect_DIR=${build_dir}"
    quit_if_fail "${FUNCNAME[0]}: cmake inside ${plugin_to_dir} failed"
    echo "make"
    eval "make"
    quit_if_fail "${FUNCNAME[0]}: make inside ${plugin_to_dir} failed"
}

fix_test(){
    # Inputs:
    #   $1: a folder of build
    #   $2: a test
    local aspect_dir="$1"
    local build_dir="$2"
    local name_of_test="$3"
    # fix all files
    test_source_dir="${aspect_dir}/tests/${name_of_test}"
    test_output_dir="${aspect_dir}/${build_dir}/tests/output-${name_of_test}"
    [[ -d ${test_source_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ${test_source_dir} doesn't exist"; exit 1; }
    [[ -d ${test_output_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ${test_output_dir} doesn't exist"; exit 1; }
    echo "Fixing test $3"
    for _file in ${test_source_dir}/*; do
        file_name=$(basename "${_file}")
        _file_out="${test_output_dir}/${file_name}"
        [[ -e ${_file_out} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ${_file_out} doesn't exist, please run the test first"; exit 1; }
        cp "${_file_out}" "${test_source_dir}"
    done
    # redo the test
    local current_dir=$(pwd)
    cd "${aspect_dir}/${build_dir}"
    echo "remake tests"
    eval "make setup_tests"
    echo "redo test ${name_of_test}"
    eval "ctest -R ${name_of_test}"
    cd "${current_dir}"
}


main(){
    ###
    # main function
    ###
    [[ -n "$1" ]] && command="$1" || cecho $BAD "\$1 must be given for options"
    if [[ "${command}" = "-h" ]]; then
        usage
    elif [[ "${command}" = "all" ]]; then
        ##
        # Build the main program with all the plugins
        # Inputs:
        # Terninal Outputs
        ##
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a name of folder"
        build_aspect_project "$2" "release"
        build_aspect_project "$2" "debug"
    elif [[ "${command}" = "all_plugins" ]]; then
        ##
        # Build all the plugins separately
        # Inputs:
        # Terninal Outputs
        ##
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a name of folder"
        build_all_plugins_separately "$2"
    elif [[ "${command}" = "fix_test" ]]; then
        ##
        # todo
        # Fix test
        # Inputs:
        # Terninal Outputs
        ##
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a test source folder (aspect/tests)"
    	[[ -n "$3" ]] || cecho $BAD "\$3 must be a name of build folder"
    	[[ -n "$4" ]] || cecho $BAD "\$4 must be a name of test"
        fix_test "$2" "$3" "$4"
    else
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a name of folder"
        local build_dir="${ASPECT_SOURCE_DIR}/build_$2"
        [[ -z $3 ]] && build_aspect_plugin "${build_dir}" "${command}" || build_aspect_plugin "${build_dir}" "${command}" "$3"
    fi
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
