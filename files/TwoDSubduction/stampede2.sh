#!/bin/bash
#SBATCH -J 2d_subduction           # Job name
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -p skx-normal         # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 48              # Total # of mpi tasks
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=hylli@ucdavis.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

module list
pwd
date


# Launch MPI code...

ibrun "${ASPECT_SOURCE_DIR}/build_TwoDSubduction/aspect" "case.prm"
