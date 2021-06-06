#!/bin/bash -l
#SBATCH -J task
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --threads-per-core=2
#SBATCH --tasks-per-node=64
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 96:00:00
#SBATCH --partition=high2

module unload openmpi
export PATH=/home/rudolph/sw/openmpi-4.0.5/bin:$PATH

srun  /home/lochy/software/aspect/build_TwoDSubduction/aspect case.prm
