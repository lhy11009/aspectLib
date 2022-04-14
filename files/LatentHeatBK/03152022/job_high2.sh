#!/bin/bash -l
#SBATCH -J task
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --threads-per-core=1
#SBATCH --tasks-per-node=8
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 24:00:00
#SBATCH --partition=high2
#SBATCH --switches=1

module unload openmpi/4.0.1
module load openmpi/4.1.0-mpi-io
>&2 echo "list of modules:"
>&2 module list
>&2 echo "aspect source: ${ASPECT_SOURCE_DIR}"

srun  ${ASPECT_SOURCE_DIR}/build_master_TwoD/aspect case.prm
