#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir}/utilities.sh"

filename="test.prm"
nnode=1
total_tasks=1
time_by_hour=24
partition="med2"
name="task"
mem_per_cpu=2000  # 2000M

test_dir="${dir}/.test"  # do test in this directory
if ! [[ -d ${test_dir} ]]; then
    mkdir "${test_dir}"
fi
test_fixtures_dir="${dir}/tests/integration/fixtures"

usage()
{
    printf "\
Submit a job to cluster with a slurm system

Usage:
  %s [options] [server_name] [file_name]

Options:
"
}

parse_command(){
    unset _command
    # parse command with the fist input
    if [[ "$1" =~ -.* || "$1" = "submit" ]]; then
	    _command="submit"
    elif [[ "$1" = "remote" ]]; then
	    _command="remote"
    elif [[ "$1" = "test" ]]; then
	    _command="test"
    elif [[ "$1" = "remote_test" || "$1" = "rt" ]]; then
	    _command="remote_test"
    elif [[ "$1" = "write_log" ]]; then
	    _command="write_log"
    else
	    return 1
    fi
    return 0
}

parse_options(){
    # parse parameters from command line
    # todo pass in name of case
    while [ -n "$1" ]; do
      param="$1"
      case $param in
        -h|--help)
          echo ""  # help information
          exit 0
        ;;
        #####################################
        # filename
        #####################################
        [^-]*)
          shift
          filename="$param"
          filename=${filename/#\.\//}
        ;;
        #####################################
        # number of total tasks
        #####################################
        -n)
          shift
          total_tasks="${1}"
        ;;
        -n=*|--total_tasks=*)
          total_tasks="${param#*=}"
        ;;
        #####################################
        # number of nodes
        #####################################
        -N)
          shift
          nnode="${1}"
        ;;
        -N=*|--nnode=*)
          nnode="${param#*=}"
        ;;
        #####################################
        # time in hour
        #####################################
        -t)
          shift
          time_by_hour="${1}"
        ;;
        -t=*|--time=*)
          time_by_hour="${param#*=}"
        ;;
        #####################################
        # partition
        #####################################
        -p)
          shift
          partition="${1}"
        ;;
        -p=*|--partition=*)
          partition="${param#*=}"
        ;;
        #####################################
        # memory per cpu
        #####################################
        -m)
          shift
          mem_per_cpu="${1}"
        ;;
        -m=*|--mem-per-cup=*)
          mem_per_cpu="${param#*=}"
        ;;
      esac
      shift
    done
}

submit(){
    # The useful part of your job goes below

    export OMP_NUM_THREADS=$SLURM_NTASKS

    # Aspect executable

    Aspect_executable="${Aspect_DIR}/aspect"
    prm_file="$dir/$filename"


    # compose the sbatch file
    # todo: add comment to sbatch file

    if [ -f 'job.sh' ]; then
    	eval "rm job.sh"
    fi
    eval "touch job.sh"
    echo "#!/bin/bash -l" >> job.sh
    echo "#SBATCH -J $name" >> job.sh
    echo "#SBATCH -N $nnode" >> job.sh
    echo "#SBATCH -n $total_tasks" >> job.sh
    echo "#SBATCH -o $name-%j.stdout" >> job.sh
    echo "#SBATCH -e $name-%j.stderr" >> job.sh
    echo "#SBATCH -t $time_by_hour:00:00" >> job.sh
    echo "#SBATCH --partition=$partition" >> job.sh
    echo "#SBATCH --mem-per-cpu=$mem_per_cpu" >> job.sh
    echo "" >> job.sh
    echo "export OMP_NUM_THREADS=\$SLURM_NTASKS" >> job.sh
    echo "" >> job.sh
    echo "srun ${Aspect_executable} ${prm_file}" >> job.sh

    # submit the job

    eval "sbatch -p $partition job.sh"
}





################################################################################
# test functions
################################################################################
test_parse_command(){
    unset _command
    # case 1
    parse_command 'submit' '-n' '32'
    if [[ ${_command} != "submit" ]]; then
	cecho ${BAD} "test_parse_command fail for \"submit -n 32\""
        exit 1
    fi
    unset _command
    # case 2
    parse_command '-n' '32'
    if [[ ${_command} != "submit" ]]; then
	cecho ${BAD} "test_parse_command fail for \"-n 32\""
	exit 1
    fi
    unset _command
    # case 3
    parse_command 'remote' '-n' '32'
    if [[ ${_command} != "remote" ]]; then
	cecho ${BAD} "test_parse_command fail for \"remote -n 32\""
	exit 1
    fi
    unset _command
    # case 4
    parse_command 'foo'
    if [[ "$?" != 1 ]]; then
	cecho ${BAD} "test_parse_command fail for \"foo\""
	exit 1
    fi
    cecho ${GOOD} "test_parse_command passed"
}

test_submit(){
    # test submit to local slurm system
    cd "${test_dir}"  # shift to testdir
    local job_id
    # test 1: mission with 1 core
    local _test="slurm.sh -n 1 -p med2 ${dir}/tests/integration/fixtures/submit_test.prm"
    job_id=$(eval "${_test}" | sed 's/Submitted\ batch\ job\ //')
    if ! [[ ${job_id} =~ ^[0-9]*$ ]]; then
        cecho ${BAD} "test_submit fail for \"${_test}\", job id is not returned"
	exit 1
    fi
    # get info from the squeue command
    get_job_info ${job_id} 'ST'
    quit_if_fail "get_job_info: no such stat 'ST'"
    if ! [[ ${return_value} = 'R' || ${return_value} = 'PD' ]]; then
        cecho ${BAD} "test_submit fail for \"${_test}\", job failed"
	exit 1
    fi
    eval "scancel ${job_id}"  # terminate this job

    # test 2: mission with nproc core
    local _nproc=$(nproc)  # numbers of cores in a node
    _test="slurm.sh -n ${_nproc}  -p med2 ${dir}/tests/integration/fixtures/submit_test.prm"
    job_id=$(eval "${_test}" | sed 's/Submitted\ batch\ job\ //')
    if ! [[ ${job_id} =~ ^[0-9]*$ ]]; then
	cecho ${BAD} "test_submit fail for \"${_test}\", job id is not returned"
	exit 1
    fi
    # get info from the squeue command
    get_job_info ${job_id} 'ST'
    quit_if_fail "get_job_info: no such stat 'ST'"
    if ! [[ ${return_value} = 'R' || ${return_value} = 'PD' ]]; then
	cecho ${BAD} "test_submit fail for \"${_test}\", job failed"
	exit 1
    fi
    get_job_info ${job_id} 'CPU'
    quit_if_fail "get_job_info: no such stat 'CPU'"
    if ! [[ ${return_value} = ${_nproc} ]]; then
	cecho ${BAD} "test_submit fail for \"${_test}\", nproc doesn't match"
	exit 1
    fi
    eval "scancel ${job_id}"  # terminate this job
    cecho ${GOOD} "test_submit passed"
    cd "${dir}"  # shift back to main directory
}


main(){
    parse_command "$1" # parse the command
    quit_if_fail "No such command \"$1\""
    if [[ ${_command} = "test" ]]; then
        test_parse_command
	test_submit
    elif [[ ${_command} = "submit" ]]; then
    	parse_options "$@"  # parse option with '-'
    	submit  # submit job
    elif [[ ${_command} = "remote_test" ]]; then
	    # test submit to cluster
        local server_info="$2"  # user@server
        shift
        if ! [[ ${server_info} =~ .*\@.* ]]; then
            # check format
            cecho ${BAD} "with 'remote' command, '\$2' needs to be 'user@server'"
        fi
        # submit job
        # todo, drag down results and compare
        # todo reset the scheme of output characters on server
        ssh ${server_info} << EOF
    	    eval "slurm.sh test"
EOF
    else
        # command error is already catched
        echo "foo"
    fi
}


if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
	main $@
fi
