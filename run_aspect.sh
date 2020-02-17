#!/bin/bash
# submitjob locally

dir='.'
name='test'
executable="${ASPECT_SOURCE_DIR}/build/aspect"
total_tasks=1

while [ -n "$1" ]; do
    param="$1"
    case $param in

        -h|--help)
            echo "Submit aspect job locally
Usege: run_aspect.sh [filename]

-n, --total_tasks=
    total parallel tasks to run" # help information
            exit 0
        ;;
        #####################################
        # filename
        #####################################
        [^-]*)
            filename="$param"
            filename=${filename/#\.\//}
        ;;
        #####################################
        # number of parallel tasks
        #####################################
        -n)
            shift
            total_tasks="${1}"
        ;;
        -n=|--total_tasks=*)
            total_tasks="${param#*=}"
        ;;
    esac
    shift # go to next variable, so we can still do something
done

cd ${dir}
mpirun -np $total_tasks $executable $filename >"job.stdout" 2>"job.stderr" &
