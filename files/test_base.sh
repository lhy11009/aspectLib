#!/bin/bash

################################################################################
# Tests functions in ().sh
# Run:
#   ./().sh
# Dependency:
#   env:
#       ASPECT_LAB_DIR
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"


test_foo(){
    ###
    # (descprition)
    ###
    # compare_files "${FUNCNAME[0]}" ${job_sh_file_std} ${job_sh_file}  # compile file contents
}


main(){
    ###
    # ()
    ###
    local local_passed_tests=0
    local local_failed_tests=0
    
    # test 1
    test_foo
    if [[ $? -eq 0 ]]; then
        ((local_passed_tests++))
    else
        ((local_failed_tests++))
    fi

    # message
    final_message "test_().sh" ${local_passed_tests} ${local_failed_tests}

    return 0
}

main $@
