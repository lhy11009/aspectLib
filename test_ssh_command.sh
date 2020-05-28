# !/bin/bash

source record.sh

FILENAME="test.prm"  # filename for aspect


##
# Print how to use

usage()
{
    printf "\
Submit a job to cluster with a slurm system

Usage:
  %s [options] [server_name] [file_name]

Options:
"
}

server_submit_job(){
    ssh $1 << EOF
        eval "cd $(dirname $2)"
        take_record 'cd $(dirname $2)' '$HOME/server_runs'
        eval "submit_job.sh $(basename $2)"
        take_record 'submit_job.sh $(basename $2)' '$HOME/server_runs'
EOF
}

main()
{
    while [ -n "$1" ]; do
        param="$1"
        case $param in
            -h|--help)
                usage
                exit 0
            ;;
            #####################################
            # Prefix path
            -p)
                shift # shift seems to move all parameters forward
                PREFIX="${1}"
            ;;
            -p=*|--prefix=*)
                PREFIX="${param#*=}" # '#' hear seems to eliminate '*=' from param, haoyuan
                # replace '~' by $HOME
                PREFIX=${PREFIX/#~\//$HOME\/}
            ;;
        esac
        shift # go to next variable, so we can still do something
    done
    server_submit_job 'lochy@peloton.cse.ucdavis.edu' '\$TwoDSubduction_DIR/katrina_case/katrina_case_parse_inputs_1/test.prm'
}

main "$@"
