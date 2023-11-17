#!/bin/bash

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
    echo "cmake -DAspect_DIR=${build_dir} ."
    eval "cmake -DAspect_DIR=${build_dir} ."
    echo "make"
    eval "make"
}


################################################################################
# run the docker
################################################################################
# paths of folders & check they exist
ASPECT_SOURCE_DIR="/home/dealii/aspect"
WORLD_BUILDER_SOURCE_DIR="/home/worldbuilder"
ASPECT_BUILD_DIR="${ASPECT_SOURCE_DIR}/tester-build"
[[ -d ${ASPECT_SOURCE_DIR} ]] || { cecho ${BAD} "${FUNCNAME[0]}: ASPECT_SOURCE_DIR doesn't exist"; exit 1; }
[[ -d ${WORLD_BUILDER_SOURCE_DIR} ]] || { cecho ${BAD} "${FUNCNAME[0]}: WORLD_BUILDER_SOURCE_DIR doesn't exist"; exit 1; }
[[ -d ${ASPECT_BUILD_DIR} ]] && rm -r ${ASPECT_BUILD_DIR}
mkdir ${ASPECT_BUILD_DIR}

# add safe to dubious ownership in repository at '/home/dealii/aspect'
git config --global --add safe.directory /home/dealii/aspect

# build the main program
cd ${ASPECT_BUILD_DIR}
cmake -G "Ninja" -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON -D WORLD_BUILDER_SOURCE_DIR="${WORLD_BUILDER_SOURCE_DIR}" ..
ninja -j 6

# setup tests
ninja setup_tests

# build plugins
build_plugin_dir="${ASPECT_SOURCE_DIR}/tester-build"
plugins_dir="${ASPECT_SOURCE_DIR}/plugins"
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
# build
for plugin in ${plugins[@]}; do
    build_aspect_plugin "${ASPECT_BUILD_DIR}" "${plugin}" "${build_plugin_dir}"
done
