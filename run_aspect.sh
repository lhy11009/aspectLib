#!/bin/bash
# submitjob locally

dir='.'
name='test'
total_tasks=1

while [ -n "$1" ]; do
    param="$1"
    case $param in

        -h|--help)
            echo "Submit aspect job locally
Usege: run_aspect.sh [filename]

-n, --total_tasks=
    total parallel tasks to run

-p, --project=
    project name, so the the executable will be:
        'build_\${project}/aspect', otherwise it is 'build/aspect'

example usages:
    run_aspect.sh -p TwoDSubduction -n 4 case.prm " # help information
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
        #####################################
        # project_folder
        #####################################
        -p)
            shift
            project="${1}"
        ;;
        -p=|--project=*)
            project="${param#*=}"
        ;;
    esac
    shift # go to next variable, so we can still do something
done

# fix executable
[[ -z ${project} ]] && executable="${ASPECT_SOURCE_DIR}/build/aspect" || executable="${ASPECT_SOURCE_DIR}/build_${project}/aspect"

cd ${dir}
mpirun -np $total_tasks $executable $filename >"job.stdout" 2>"job.stderr" &
