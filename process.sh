#!/bin/bash

source server_commands.sh

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
        if [[ ! ${#ARGUMENT[@]} -eq 1 ]]; then
            exit "With the submit option, Argumemnt mush be just the name for the file on server"
        fi
        REMOTE_FILE=${ARGUMENT[0]}
        # REMOTE_FILE='\$TwoDSubduction_DIR/katrina_case/katrina_case_parse_inputs_1/test.prm' # test remote file
        server_submit_job "$USER@$SERVER" $REMOTE_FILE
    fi
}

main "$@"
