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
    # The useful part of your job goes below

    export OMP_NUM_THREADS=$SLURM_NTASKS

    # Aspect executable
    # add project_dir after aspect_dir, like 'build_master'

    Aspect_project_DIR="${Aspect_DIR}_${project}"
    Aspect_executable="${Aspect_project_DIR}/aspect"

    # Aspect_executable="${Aspect_DIR}/aspect"

    # compose the sbatch file
    # future: add comment to sbatch file

    # check for filename

    [[ -e "${filename}" ]] || { cecho ${BAD} "${filename} doesn't exit"; exit 1; }

    # first cd to the folder
    local previous_dir=$(pwd)
    local case_dir=$(dirname "$filename")
    cd "$case_dir"

    if [ -f 'job.sh' ]; then
    	eval "rm job.sh"
    fi

    # create slurm script
    eval "touch job.sh"

    # write slurm messages
    echo "#!/bin/bash -l" >> job.sh
    echo "#SBATCH -J $name" >> job.sh
    echo "#SBATCH -N $nnode" >> job.sh
    echo "#SBATCH -n $total_tasks" >> job.sh
    # tasks per node
    ((tasks_per_nodes=total_tasks/nnode))
    echo "#SBATCH --threads-per-core=2" >> job.sh
    echo "#SBATCH --tasks-per-node=${tasks_per_nodes}" >> job.sh
    echo "#SBATCH -o $name-%j.stdout" >> job.sh
    echo "#SBATCH -e $name-%j.stderr" >> job.sh
    echo "#SBATCH -t $time_by_hour:00:00" >> job.sh
    echo "#SBATCH --partition=$partition" >> job.sh
    echo "#SBATCH --mem-per-cpu=$mem_per_cpu" >> job.sh
    echo "" >> job.sh

    # unload module openmpi and load Max Rudolph's version
    if [[ $(hostname) = "peloton.cse.ucdavis.edu" ]]; then
      echo "module unload openmpi" >> job.sh
      echo "export PATH=/home/rudolph/sw/openmpi-4.0.5/bin:\$PATH" >> job.sh
      echo "" >> job.sh
    fi

    addition=""
    [[ -n ${bind_to} ]] && addition="$addition --cpu-bind=${bind_to}"
    echo "srun ${addition} ${Aspect_executable} ${filename}" >> job.sh

    # submit the job, hold if the hold option is 1

    (( hold == 0 )) && eval "sbatch -p $partition job.sh"

    # go back to previous dir
    cd "${previous_dir}"
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
