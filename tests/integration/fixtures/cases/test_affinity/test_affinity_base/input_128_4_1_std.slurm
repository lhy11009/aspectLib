#!/bin/bash
#SBATCH -p rome-256-512
#SBATCH -n 128
#SBATCH --exclusive
#SBATCH --ntasks-per-node=128
#SBATCH --time=02:00:00
#SBATCH --job-name=run_128_4_1
#SBATCH --switches=1
set -x
module unload openmpi/4.0.1
            module load openmpi/4.1.0-mpi-io
module list
ulimit -l unlimited
mpirun -n 128 --bind-to socket --report-bindings /home/lochy/Softwares/aspect/build_master_TwoD/aspect .test/test_affinity/test_affinity_base/tmp/peloton-rome-128tasks-socket-openmpi-4.1.0/input_128_4_1
