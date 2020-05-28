#!/bin/bash
# init new case on cluster by uploading a .prm file

PROJECT_DIR="/home/lochy/ASPECT_PROJECT"
PROJECT_DIR_ON_CLUSTER="/home/lochy/ASPECT_PROJECT"

init_new_case_on_cluster(){
    FOLDER=$1
    CASE_NAME=$2
    CLUSTER=$3
    COMMAND="rsync -avur $PROJECT_DIR/$FOLDER/$CASE_NAME/*.prm $CLUSTER:$PROJECT_DIR_ON_CLUSTER/$FOLDER/$CASE_NAME/"
    echo "command: $COMMAND"
    eval $COMMAND
}

init_new_case_on_cluster "DATA/katrina_case" "katrina_case_parse_inputs_4" "lochy@peloton.cse.ucdavis.edu"

# eval "python -m shilofue -p twoDSubduction -t parse -i files/twoDSubduction_base.prm -o test.prm"
