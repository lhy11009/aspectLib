#!/bin/bash
#look at jobs in mpi
only_main=0
while [ -n "$1" ]; do
    param="$1"
    case $param in

        -h|--help)
            echo "show aspect jobs in system
-m, --main
    only show main process"  # help information
            exit 0
        ;;
        #####################################
        # Prefix path
        -m|--main)
            shift # shift seems to move all parameters forward
            only_main=1
        ;;
    esac
    shift # go to next variable, so we can still do something
done

if [ $only_main -eq 0 ]; then
	ps -aux | grep build.*aspect # show all process
else
	ps -aux | grep mpirun.*build.*aspect # only show main process
fi
