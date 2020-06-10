#!/bin/bash -l

# Name of the job
#SBATCH -J test

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o %j.stdout
#SBATCH -e %j.stderr

# envirmont variables below

dir=$(pwd)
filename="test.prm"
nnode=1
total_tasks=1
time_by_hour=24
partition="med2"
name="task"
mem_per_cpu=2000  # 2000M

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
    exit 0
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

main(){
    parse_options  # parse options with '-'
    submit  # submit job
}
