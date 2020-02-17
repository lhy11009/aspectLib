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
Nproc_per_node=1

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
  esac
  shift
done 

echo $filename
# The useful part of your job goes below

export OMP_NUM_THREADS=$SLURM_NTASKS

# Aspect executable

Aspect_executable="${Aspect_DIR}/aspect"
prm_file="$dir/$filename"

# The main job executable to run: note the use of srun before it
srun -n $total_tasks -N $nnode ${Aspect_executable} ${prm_file}
