#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --threads-per-core=1
#SBATCH --tasks-per-node=128
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 48:00:00
#SBATCH -A billen
#SBATCH --partition=p-billen
#SBATCH --mem-per-cpu=4G
#SBATCH --switches=1

module unload openmpi
module load gcc/11.3.0
module load openmpi/4.1.5-mpi-io
module load cmake/3.26.3

source /group/billengrp/Software/deal.ii/deal.ii-9.5.0-Native-32bit-candi-gcc-11.3.0-openmpi4.1.5-mpi-io-rome-256-512/configuration/enable.sh

>&2 echo "list of modules:"
>&2 module list
>&2 echo "aspect source: ${ASPECT_SOURCE_DIR}"

srun  ${ASPECT_SOURCE_DIR}/build_master_TwoD_p-billen/aspect case.prm
