#!/bin/bash -l
#SBATCH -J task
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --tasks-per-node=16
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 24:00:00
#SBATCH --partition=high2
#SBATCH --mem-per-cpu=2000

module unload openmpi
export PATH=/home/rudolph/sw/openmpi-4.0.5/bin:$PATH

srun  /home/lochy/software/aspect/build_TwoDSubduction/aspect /home/lochy/ASPECT_PROJECT/TwoDSubduction/slurm_test_submit2/case.prm
