# !/bin/bash

FILENAME="test.prm"  # filename for aspect
ssh lochy@peloton.cse.ucdavis.edu << EOF
    eval "cd \$TwoDSubduction_DIR/katrina_case/katrina_case_parse_inputs_1"
    eval "pwd > ~/test_results"
    if [[ -n "$FILENAME" && -e "$FILENAME" ]]; then
        eval "submit_job.sh $FILENAME"
    fi
EOF
