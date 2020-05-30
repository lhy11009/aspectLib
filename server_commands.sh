# !/bin/bash

##
# Print how to use


server_submit_job(){
    ssh $1 << EOF
        eval 'source \$ASPECT_LAB_DIR/record.sh'
        eval "cd $(dirname $2)"
        take_record 'cd $(dirname $2)' '\$HOME/server_runs'
        eval "submit_job.sh $(basename $2)"
        take_record 'submit_job.sh $(basename $2)' '\$HOME/server_runs'
EOF
}

export -f server_submit_job
