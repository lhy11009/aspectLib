#!/bin/bash
# case manager
# Usage:
#   ./aspect_lib.sh + command + options

main(){
    # parameter list, todo
    # case_name =
    if [[ ${_commend} = 'create' ]]; then
        python -m shilofue.TwoDSubduction create $case_name
        # assertion todo
        # assert()
        # submit to server
        # slurm.sh ..
    elif [[ ${_commend} = 'create_group' ]]; then
        python -m shilofue.TwoDSubduction create_group ${group_name}
        # assertion todo
        # assert()
        # submit to server
        # slurm.sh ..
    elif [[ ${_commend} = 'test' ]]; then
        # do test, todo
        echo "foo"
    fi
    return 0
}

main $@
