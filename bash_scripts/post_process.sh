#!/bin/bash


################################################################################
# Post process on cases (TODO)
#
# Dependencies:
#    env:
#        ASPECT_LAB_DIR
#        ASPECT_SOURCE_DIR
#    scripts:
#        This script depends on other post-process scripts(see below):
#
# Example Usage:
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
shilofue_dir = "${ASPECT_LAB_DIR}/shilofue"
py_run_time = "${shilofue_dir}/PlotRunTime.py"
py_depth_average = "${shilofue_dir}/PlotDepthAverage.py"
# source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"


#parse_options(){
#    ###
#    # parse options
#    ###
#    while [ -n "$1" ]; do
#      param="$1"
#      case $param in
#        -h|--help)
#          usage  # help information
#          exit 0
#        ;;
#      esac
#      shift
#    done
#}

TwoDSubduction_post_process(){
    ###
    # post process for TwoDSubduction
    # Inputs:
    #   $1: case_dir
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    ###
    local case_dir = "$1"
    # go to parental directory
    local _dir = $(pwd)
    cd "${ASPECT_LAB_DIR}"
    # run time
    eval "python -m shilofue.PlotRunTime plot_case -i ${case_dir}"
    # newton solver
    eval "python -m shilofue.PlotRunTime"
    # depth average
    local depth_average_path = "${case_dir}/output/depth_average.txt"
    eval "python -m shilofue.PlotDepthAverage plot_by_time -i ${depth_average_path}"

    # go back to original directory
    cd "${_dir}"

    return 0
}


main(){
    ###
    # main function
    ###
    if [[ "$1" = "TwoDSubduction_case" ]]; then
        ##
        # (Descriptions)
        # Innputs:
        # Terninal Outputs
        ##
        shift
        case_dir = "$1"
        TwoDSubduction_post_process "${case_dir}"
    else
	    cecho "${BAD}" "option ${1} is not valid\n"
    fi
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi

## notes

#trouble shooting
# [[ -n "$2" ]] || { cecho "${BAD}" "no log file given ($2)"; exit 1; }

#debuging output
# printf "${FUNCNAME[0]}, return_values: ${return_values[@]}\n" # debug

# parse options
        # shift
        # parse_options() $@
