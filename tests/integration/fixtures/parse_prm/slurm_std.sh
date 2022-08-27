##!/bin/bash -l
#SBATCH -J task
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --threads-per-core=2
#SBATCH --tasks-per-node=64
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 300:00:00
#SBATCH -A billen
#SBATCH --partition=p-billen
#SBATCH --switches=1
module unload openmpi/4.0.1
module unload deal.II
module load openmpi/4.1.0-mpi-io
module load /group/billengrp/software/deal.ii/deal.ii-9.3.3-Native-32bit-candi-gcc-11.1.0-openmpi4.1.0-mpi-io-rome-256-512/configuration/modulefiles/default
srun ${ASPECT_SOURCE_DIR}/build_master_TwoD/aspect case.prm