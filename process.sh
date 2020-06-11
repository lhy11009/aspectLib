#!/bin/bash

# todo
# complete usage
usage()
{
    printf "\

Usage:
Run command on the server side

    %s [server_name] [command] [file_name]

server_name:

    Name of server

command:

    submit:
        Submit a job to cluster with a slurm system
"
}

# Submit job to servers
# Inputs:
#   $1(str): user@server
#   $2(str): filename for .prm
#   $3(str): addtional variables
process_submit(){
    get_addtional_options $3  # get addtional output
    ssh $1 << EOF
        eval 'source \$ASPECT_LAB_DIR/utilities.sh'
        eval "cd $(dirname $2)"
        take_record 'cd $(dirname $2)' '\$HOME/server_runs'
        eval "submit_job.sh $ADDITIONAL_OPTIONS $(basename $2)"
        take_record 'submit_job.sh $ADDITIONAL_OPTIONS $(basename $2)' '\$HOME/server_runs'
EOF
    unset ADDITIONAL_OPTIONS
}

# get_addtional_options from a string which looks like:
# '#--total_tasks=32#--partition=med2'
# and return with something like:
# ' --total_tasks=32 --partition=med2'
get_addtional_options(){
    ADDITIONAL_OPTIONS=${1//\#/ }
}

# Main function
main()
{
    unset ARGUMNET
    unset USER
    unset SERVER
    unset OPTION
    unset CASE_NAME
    unset DIR
    while [ -n "$1" ]; do
        param="$1"
        case $param in
            -h|--help)
                usage
                exit 0
            ;;
            #####################################
            # Arguments
            [^-]*)
                ARGUMENT="${ARGUMENT} ${1}"
            ;;
            #####################################
            # Option
            -o)
                shift
                OPTION="${1}"
            ;;
            -o=*|--option=*)
                OPTION="${param#*=}"
            ;;
            #####################################
            # Casename
            -c)
                shift
                CASE_NAME="${1}"
            ;;
            -c=*|--case_name=*)
                CASE_NAME="${param#*=}"
            ;;
            #####################################
            # Path to case data
            -d)
                shift
                DIR="${1}"
            ;;
            -d=*|--dir=*)
                DIR="${param#*=}"
                # replace '~' by $HOME
                DIR=${PREFIX/#~\//$HOME\/}
            ;;
            #####################################
            # USER
            -u)
                shift
                USER="${1}"
            ;;
            -u=*|--user=*)
                USER="${param#*=}"
            ;;
            #####################################
            # Server name
            -s)
                shift
                SERVER="${1}"
            ;;
            -s=*|--server=*)
                SERVER="${param#*=}"
            ;;
            #####################################
        esac
        shift # go to next variable, so we can still do something
    done
    # handle cases with wrong inputs
    if [[ ! -v "OPTION" || ! -n "OPTION" ]]; then
        exit "Must give an OPTION variable"
    fi
    # REMOTE_DIR=  # get the remote path, todo
    if [[ "$OPTION" = "submit" ]]; then
        # submit a case to remote server
        if [[ ! -v "USER" || ! -n "USER" ]]; then
            exit "With the submit option, Must give a USER variable"
        fi
        if [[ ! -v "SERVER" || ! -n "SERVER" ]]; then
            exit "With the submit option, Must give a SERVER variable"
        fi
        if [[ ${#ARGUMENT[@]} -gt 2 ]]; then
            exit "With the submit option, Argumemnt mush be just the name for the file on server, an addtional one may be given as additional option"
        fi
        local REMOTE_FILE=${ARGUMENT[0]}
        local ADDITIONAL_OPTIONS=${ARGUMENT[1]}
        # REMOTE_FILE='\$TwoDSubduction_DIR/katrina_case/katrina_case_parse_inputs_1/test.prm' # test remote file
        process_submit "$USER@$SERVER" $REMOTE_FILE $ADDITIONAL_OPTION
    fi
}

main "$@"
