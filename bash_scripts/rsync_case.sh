#!/bin/bash


################################################################################
# rsync from server
#
# Dependencies:
#    env:
#        ASPECT_LAB_DIR
#        ASPECT_SOURCE_DIR
#
# Example Usage:
#    sync from remote:
#        ./bash_scripts/rsync_case.sh peloton TwoDSubduction non_linear30/eba_re
#    sync to remote:
#        ./bash_scripts/rsync_case.sh peloton TwoDSubduction non_linear30/eba_re -f true
################################################################################

dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null 2>&1 && pwd  )"
source "${ASPECT_LAB_DIR}/bash_scripts/utilities.sh"
SERVER_LIST="${ASPECT_LAB_DIR}/env/server_list"  # in this list, one need to include the address of servers they are using

usage(){
  # TODO
  # usage of this script
    _text="
${BASH_SOURCE[0]}

Dependencies:
   env:
       ASPECT_LAB_DIR
       ASPECT_SOURCE_DIR

Example Usage:
   sync from remote (add skip_restart would skip restart files):
       ./bash_scripts/rsync_case.sh peloton TwoDSubduction non_linear30/eba_re --skip_restart 1
   sync to remote:
       ./bash_scripts/rsync_case.sh peloton TwoDSubduction non_linear30/eba_re -f true

   sync aspect work place
       ./bash_scripts/rsync_case.sh ucd TwoDSubduction
"
    printf "${_text}"

}

parse_options(){
    # parse options
    while [ -n "$1" ]; do
      param="$1"
      case $param in
        -h|--help)
          usage  # help information
          exit 0
        ;;
        -f|--from_local)
          shift
          g_from_local="$1"
        ;;
        --skip_restart)
          shift
          skip_restart_file="$1"
        ;;
      esac
      shift
    done
}

rsync_server(){
    ###
    # rsync from peloton
    # Inputs:
    # Global variables as input:
    # Return:
    # Global variables as output(changed):
    ###
    local case_dir=$1
    local local_dir=$2
    local server_dir=$3
    local addr=$4
    local flag=""
    [[ -n "${case_dir}" ]] || { cecho "${BAD}" "no case given (\$1}"; exit 1; }
    [[ -n "${local_dir}" ]] || { cecho "${BAD}" "fail to read local directory"; exit 1; }
    [[ -n "${server_dir}" ]] || { cecho "${BAD}" "fail to read peloton directory"; exit 1; }
    [[ -n "${skip_restart_file}" ]] && flag="--exclude=*restart*" 
    if [[ -n "${g_from_local}" ]]; then
        local parental_dir=$(dirname "${server_dir}/${case_dir}")
        echo "rsync -avur --progress ${flag} ${local_dir}/${case_dir} ${addr}:${parental_dir}/"
        eval "rsync -avur --progress ${flag} ${local_dir}/${case_dir} ${addr}:${parental_dir}/"
    else
        local parental_dir=$(dirname "${local_dir}/${case_dir}")
        echo "rsync -avur --progress ${flag} ${addr}:${server_dir}/${case_dir} ${parental_dir}/"
        eval "rsync -avur --progress ${flag} ${addr}:${server_dir}/${case_dir} ${parental_dir}/"
    fi
    return 0
}

rsync_workplace(){
    ### 
    # rsync the work place from server(figures, docs)
    ###
    local root_dir=$1
    [[ -n "${root_dir}" ]] || { cecho "${BAD}" "no root directory given (\$1}"; exit 1; }

    # project, e.g. TwoDSubduction
    local project=$3

    echo "rsync -avur --progress $2/* --include=\"case.prm\" --include=\"log.txt\" --include=\"*img/*\" --include=\"*/\" --exclude=* ${root_dir}/"
    eval "rsync -avur --progress $2/* --include=\"case.prm\" --include=\"log.txt\" --include=\"*img/*\" --include=\"*/\" --exclude=* ${root_dir}/"
    echo "rsync -avur --progress ${root_dir}/* --include=\"case.prm\" --include=\"log.txt\" --include=\"*img/*\" --include=\"*/\" --exclude=* ${2}/"
    eval "rsync -avur --progress ${root_dir}/* --include=\"case.prm\" --include=\"log.txt\" --include=\"*img/*\" --include=\"*/\" --exclude=* ${2}/"
    if [[ "${project}" = "TwoDSubduction" ]]; then
        echo "rsync -avur --progress $2/* --include=\"*.md\" --include=\"*.mkd\" --include=\"*.json\" --include=\"particle.dat\" --include=\"*.sh\" --include=\"*/\" --exclude=* ${root_dir}/"  # mkd, particle and job
        eval "rsync -avur --progress $2/* --include=\"*.md\" --include=\"*.mkd\" --include=\"*.json\" --include=\"particle.dat\" --include=\"*.sh\" --include=\"*/\" --exclude=* ${root_dir}/"  # mkd, particle and job
    fi
}

main(){
    ###
    # main function
    ###
    if [[ "$1" = "-h" ]]; then
        usage
    elif [[ "$1" = "peloton" ]]; then
        ##
        # tranfer data to & from peloton
        # Innputs:
        # Terninal Outputs
        ##
        shift;
        project="$1" 
        shift;
        case_dir="$1"; shift
        parse_options $@
        local query="${project}_DIR"
        local_dir=${!query}
        [[ -n $local_dir ]] || { cecho "${BAD}" "Query of local_dir failed"; exit 1; }
        echo "local_dir: ${local_dir}"
        peloton_dir=`awk '{if ($2 ~ /'"${project}_DIR"'/){split($2, array, "="); gsub("\"", "", array[2]) ;print array[2]}}' "${ASPECT_LAB_DIR}/env/enable_peloton.sh"`
        [[ -n $peloton_dir ]] || { cecho "${BAD}" "Query of peloton_dir failed"; exit 1; }
        echo "peloton_dir: ${peloton_dir}"  # debug
        local peloton_addr=$(eval "awk '{ if(\$1 == \"peloton\") print \$2;}' $SERVER_LIST")
        rsync_server "${case_dir}" "${local_dir}" "${peloton_dir}" "${peloton_addr}"
    elif [[ "$1" = "stampede2" ]]; then
        ##
        # tranfer data to & from stampede2
        # Innputs:
        # Terninal Outputs
        ##
        shift;
        project="$1" 
        shift;
        case_dir="$1"; shift
        parse_options $@
        local query="${project}_DIR"
        local_dir=${!query}
        stampede2_dir=`awk '{if ($2 ~ /'"${project}_DIR"'/){split($2, array, "="); gsub("\"", "", array[2]) ;print array[2]}}' "${ASPECT_LAB_DIR}/env/enable_stampede2.sh"`
        echo "stampede2_dir: ${stampede2_dir}"  # debug
        local stampede2_addr=$(eval "awk '{ if(\$1 == \"stampede2\") print \$2;}' $SERVER_LIST")
        echo "stampede2_addr: ${stampede2_addr}"  # debug
        rsync_server "${case_dir}" "${local_dir}" "${stampede2_dir}" "${stampede2_addr}"
    elif [[ "$1" = "frontera" ]]; then
        ##
        # tranfer data to & from frontera
        # Innputs:
        # Terninal Outputs
        ##
        shift;
        project="$1" 
        shift;
        case_dir="$1"; shift
        parse_options $@
        local query="${project}_DIR"
        local_dir=${!query}
        frontera_dir=`awk '{if ($2 ~ /'"${project}_DIR"'/){split($2, array, "="); gsub("\"", "", array[2]) ;print array[2]}}' "${ASPECT_LAB_DIR}/env/enable_frontera.sh"`
        echo "frontera_dir: ${frontera_dir}"  # debug
        local frontera_addr=$(eval "awk '{ if(\$1 == \"frontera\") print \$2;}' $SERVER_LIST")
        echo "frontera_addr: ${frontera_addr}"  # debug
        rsync_server "${case_dir}" "${local_dir}" "${frontera_dir}" "${frontera_addr}"
    elif [[ "$1" = "ucd" && "$2" = "TwoDSubduction" ]]; then
        # TODO
        # rsync workplace
        local ucd_addr=$(eval "awk '{ if(\$1 == \"${1}\") print \$2;}' $SERVER_LIST")
        rsync_workplace "${TwoDSubduction_DIR}" "${ucd_addr}:${TwoDSubduction_ucd_DIR}" "$2"
    else
        ## TODO
        parse_options
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

