#!/bin/bash


################################################################################
# Build aspect branch
#
# Dependencies:
#    utilities.sh
#    env:
#        ASPECT_SOURCE_DIR
# Example Usage:
#   build main program and plugins:
#    ./build_aspect.sh all TwoDSubduction release
#   build a plugin:
#    ./build_aspect.sh subduction_temperature2d TwoDSubduction
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${dir}/utilities.sh"

build_aspect_project(){
    ###
    # build a project in aspect
    # usage of this is to bind up source code and plugins
    # Inputs:
    #   project: name of the project
    ###
    local build_dir="${ASPECT_SOURCE_DIR}/build_$1"
    [[ -d ${build_dir} ]] || mkdir ${build_dir}
    local mode
    if [[ -n $2 ]]; then
        [[ mode="debug" || mode="release" ]] || { cecho ${BAD} "${FUNCNAME[0]}: mode is either \'debug\' or \'release\'"; exit 1; }
        mode="$2"
    else
        mode="debug"
    fi

    # get the list of plugins
    plugins=("prescribe_field" "subduction_temperature2d" "slab2d_statistics" "subduction_temperature2d_ellipse")

    # build
    local current_dir=$(pwd)
    # Here we pick nproc - 1, this make sure that we don't use up all resources. 
    # But this will cause problem when nproc = 1
    # local nproc=$(($(nproc)-1))
    local nproc=8
    cd ${build_dir}
    # build source code
    eval "cmake .."
    quit_if_fail "${FUNCNAME[0]}: cmake inside ${build_dir} failed"
    eval "make ${mode}"
    quit_if_fail "${FUNCNAME[0]}: \"make ${mode}\" inside ${build_dir} failed"
    eval "make -j ${nproc}"
    quit_if_fail "${FUNCNAME[0]}: make inside ${build_dir} failed"

    # build plugins
    plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    for plugin in ${plugins[@]}; do
        build_aspect_plugin "${build_dir}" "${plugin}"
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
    
    # copy plugins
    plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
    # check plugin existence
    plugin_dir="${plugins_dir}/${plugin}"
    [[ -d ${plugin_dir} ]] || { cecho ${BAD} "${FUNCNAME[0]}: plugin(i.e. ${plugin_dir}) doesn't exist"; exit 1; }
    
    # remove old ones
    plugin_to_dir="${build_dir}/${plugin}"
    [[ -d ${plugin_to_dir} ]] && rm -r ${plugin_to_dir}
    # copy new ones
    eval "cp -r ${plugin_dir} ${build_dir}/"
    cecho ${GOOD} "${FUNCNAME[0]}: copyied plugin(i.e. ${plugin})"

    # build 
    cd ${plugin_to_dir}
    # remove cache before compling
    [[ -e "${plugin_to_dir}/CMakeCache.txt" ]] && eval "rm ${plugin_to_dir}/CMakeCache.txt"
    eval "cmake -DAspect_DIR=${build_dir}"
    quit_if_fail "${FUNCNAME[0]}: cmake inside ${plugin_to_dir} failed"
    eval "make"
    quit_if_fail "${FUNCNAME[0]}: make inside ${plugin_to_dir} failed"
}


main(){
    ###
    # main function
    ###
    [[ -n "$1" ]] && command="$1" || cecho $BAD "\$1 must be given for options"
    if [[ "${command}" = "all" ]]; then
        ##
        # Build the main program with all the plugins
        # Inputs:
        # Terninal Outputs
        ##
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a name of folder"
        build_aspect_project "$2" "$3"
    else
    	[[ -n "$2" ]] || cecho $BAD "\$2 must be a name of folder"
        local build_dir="${ASPECT_SOURCE_DIR}/build_$2"
        build_aspect_plugin "${build_dir}" "${command}"
    fi
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
