#!/bin/bash


################################################################################
# configure a slurm file
#
# Dependencies:
#    utilities.sh
#
# Example Usage:
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${dir}/utilities.sh"

configure_slurm_file(){
    ###
    # configure a slurm file
    # Inputs:
    #   $1: base file
    #   $2: target
    #   $3: options
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    ###
    [[ -n "$1" && -n "$2" ]] || cecho "${BAD}" "${FUNC_NAME[0]}: \$1 and \$2 must be given"
    local base_file="$1"
    local target_file="$2"
    # copy file
    cp "${target_file}" "${base_file}"
    # read in key and values from the file of options
    read_keys_values "$3"
    # submit key 
    i=0
    while ((i<${#keys[?]}; do
        key="${keys[i]}"
        value="${values[i]}"
        eval "sed -i \'s/${key}/${value}/g\' ${target_file}"
        ((i++))
    done
    return 0
}


main(){
    ###
    # main function
    ###
    [[ -n "$1" && -n "$2" ]] || cecho "${BAD}" "${FUNC_NAME[0]}: \$1 and \$2 must be given"
    configure_slurm_file "$1" "$2"
}

set +a  # return to default setting

##
# if this is call upon in terminal, the main function is executed
##
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
