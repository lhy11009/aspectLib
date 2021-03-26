#!/bin/bash


################################################################################
# generate a summary of a case
#
# Dependencies:
#    env:
#        ASPECT_LAB_DIR
#        ASPECT_SOURCE_DIR
#
# Example Usage:
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"

case_update(){
    ###
    # (Descriptions)
    # Inputs:
    #   $1: case directory(contains case.prm)
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    ###
    local case_dir="$1"
    local case_name=$(dirname "${case_dir}")
    prm_file="${case_dir}/case.prm"
    [[ -e "${prm_file}" ]] || { cecho ${BAD} "${FUNC_NAME[0]}: ${prm_file} doesn't exist"; return 1; }
    echo "case directory: ${case_dir}"
    # plot statistics
    echo "ploting statistics"
    # plot depth-average
    echo "ploting depth-average"
    # plot newton solver in the last 20 steps
    echo "ploting newton solver in the last 20 steps"
    # visualize last step
    echo "visualizing last step"
    # generate a markdown output
    echo "generating markdown"
    return 0
}


main(){
    ###
    # main function
    ###
    local case_dir="$1"
    case_update "${case_dir}"
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
