#!/bin/bash

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

cd ${ASPECT_BUILD_DIR}
cmake -G "Ninja" -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON -D WORLD_BUILDER_SOURCE_DIR="${WORLD_BUILDER_SOURCE_DIR}" ..

