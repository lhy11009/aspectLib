#!/bin/bash

################################################################################
# Taking record of commands
#
# 1) Example Command:
################################################################################

take_record() {
    message=$(eval echo $1)
    record_file=$(eval echo $2)

    # add time stamp to record
    timestamp () {
                echo "$(date +"%Y-%m-%d_%H-%M-%S")"

    }

    echo "$(timestamp): ${message}" >> ${record_file}

}

export -f take_record
