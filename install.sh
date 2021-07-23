#!/bin/bash


################################################################################
# Install aspectLib
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"
# subdirs
unset dir
dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
bash_subdir="${dir}/bash_scripts"
[[ -d "${bash_subdir}" ]] || { cecho "${BAD}" "install.sh: no directory (${bash_subdir})"; exit 1; }
bin_subdir="${dir}/bin"
[[ -d "${bin_subdir}" ]] || { echo "create ${bin_subdir}"; mkdir "${bin_subdir}"; }
# global variables
prefix="Lib"  # prefix to installed scripts
bashrc_outputs="" # tests to include in the .bashrc file

usage(){
  # usage of this script
    _text="
${BASH_SOURCE[0]}

Install aspectLib, run this in the containing folder!!

Dependencies:
   source files in the subdirectory of bash_scripts

Example Usage:
    installation:
        ./install.sh execute
    clean:
        ./install.sh clean
"
    printf "${_text}"

}

options(){
    ###
    # options when installing files
    ###
    list_of_bash_scripts_to_install=("rsync_case.sh" "run_aspect.sh")
}


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

install(){
    ###
    # install apsectLib
    ###
    # todo
    local dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
    options
    install_sh
    printf "${bashrc_outputs}" >> "${dir}/enable.sh" # output for bashrc
    # screen output
    printf "install.sh: installation is completed\n\
Next step: append one-liner to .bashrc:\n\
    source ${dir}/enable.sh\n"
}

clean(){
    ###
    # clean previous installation
    ###
    # todo
    # remove install executables
    printf "remove install executables\n"
    eval "[[ -d ${bin_subdir} ]] && rm ${bin_subdir}/*"
    printf "remove previous enable.sh\n"
    eval "[[ -e \"${dir}/enable.sh\" ]] && rm ${dir}/enable.sh"
}

install_sh(){
    ###
    # install bash scripts
    # Inputs:
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    #   bashrc_outputs: tests to include in the .bashrc file
    ###
    local bash_subdir="${dir}/bash_scripts"
    for bash_script_to_install in ${list_of_bash_scripts_to_install[@]}; do
        _path="${bash_subdir}/${bash_script_to_install}"
        [[ -e "${_path}" ]] || { cecho "${BAD}" "${FUNCNAME[0]} no bash scripts (${_path}), check the source files and the file list given"; exit 1; }
        _name=`echo "${bash_script_to_install}" | cut -d'.' -f1`  # todo
        _name="${prefix}_${_name}"
        _path_to="${bin_subdir}/${_name}"
        # copy
        cp "${_path}" "${_path_to}"
        # change mode
        eval "chmod +x ${_path_to}"
        cecho ${GOOD} "${FUNCNAME[0]} ${_path_to} installed"  # screen outputs
    done
    bashrc_outputs="${bashrc_outputs}\n# aspectLib executables\nexport PATH=\${PATH}:${bin_subdir}"
    return 0
}


main(){
    ###
    # main function
    ###
    if [[ "$1" = "-h" ]]; then
        usage
    elif [[ "$1" = "execute" ]]; then
        ##
        # (Descriptions)
        # Innputs:
        # Terninal Outputs
        ##
        # todo, create alias
        install
    elif [[ "$1" = "clean" ]]; then
        # clean previous installation
        clean
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
