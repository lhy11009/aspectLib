#!/bin/bash

KEY="INITIAL_ADAPTIVE_REFINEMENT"
VALUE=(1 2 4)

for VALUE in ${VALUES[@]}; do
    # todo add a python command
    # python -m shilofue.TwoDSubduction
    CASE="katrina_case_parse_inputs_1_test_$KEY_$VALUE"
    LOCAL_RPM="/home/lochy/ASPECT_PROJECT/DATA/katrina_case/$CASE/test.prm"
    REMOTE_PRM="/home/lochy/ASPECT_PROJECT/DATA/katrina_case/$CASE/test.prm"
    rsync1 $LOCAL_PRM "$user@$server:$REMOTE_PRM"
    process.sh -o submit -u lochy -s peloton.cse.ucdavis.edu "$REMOTE_PRM" ''#--total_tasks=32#--partition=med2'
done

