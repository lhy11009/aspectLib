#!/bin/bash


################################################################################
# (Descriptions)
#
# Dependencies:
#    env:
#        ASPECT_LAB_DIR
#        ASPECT_SOURCE_DIR
#
# Example Usage:
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
# source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"

some_func(){
    ###
    # (Descriptions)
    # Inputs:
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    ###
    return 0
}


main(){
    ###
    # main function
    ###
    if [[ "$1" = "foo" ]]; then
        ##
        # (Descriptions)
        # Innputs:
        # Terninal Outputs
        ##
        some_funct()
        printf "${return_values}"
    fi
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
