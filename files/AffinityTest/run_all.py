# This python script generates submission scripts (slurm-style)
# For the aspect performance benchmarks
# haoyuan: this scripts test for different number of cores as well as different level of refinement
# a tasks_per_node option is chosen, default is 32

import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt  # haoyuan: I don't think this is used

# from subprocess import run
import os

# Haoyuan: in this file, the field to change are marked with capital letters.
# e.g.  set Output directory                       = OUTPUT_DIRECTORY
base_input = "setups/spherical_shell_expensive_solver.prm"     # The 'base' input file that gets modified
#cluster_label = "PI4CS_aspect-2.0-pre-40tasks"
#cluster_label = "peloton-ii-64tasks-hwthread-openmpi-4.0.1"
#cluster_label = "peloton-ii-32tasks-core-openmpi-4.0.1"
cluster_label = "peloton-ii-32tasks-core-openmpi-4.0.1" # ?

# modify this to contain the commands necessary to setup MPI environment
#environment_setup_commands = "module load openmpi/3.1.3 intel-mkl"
# haoyuan: I think we need to load 4.0.5 here
environment_setup_commands=""


core_counts = [1,2,4,8,16,32,64,128,256,512,768,1024]#,200,300,400]#,500,800,1000,1500]
refinement_levels = [2,3,4,5]#,6]
#                                          0   1   2   3       4     5    6
minimum_core_count_for_refinement_level = [0,  0,   1,   1,   10, 100, 500]# for refinement levels 0-6
maximum_core_count_for_refinement_level = [0,  0,1000,1000, 1000,2000,2000]

setups = [1,]
tasks_per_node = 32
# make directories for temporary files
os.system('mkdir tmp')
os.system('mkdir tmp/'+cluster_label)

def generate_input_file(base_file_name,output_file_name,dictionary):
    """Read the 'base' input file from base_file_name, replace strings 
    using dictionary, and write new output file to output_file_name"""
    fh = open(base_file_name,'r')
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)  # haoyuan: run module isn't imported correctly
    # run(['rm','-f',output_file_name])
    ofh = open(output_file_name,'w')
    for line in fh:        
        for key in dictionary:
            if key in line:                
                line = line.replace(key,str(dictionary[key]))
        ofh.write(line)
    fh.close()
    ofh.close()

def generate_slurm_file(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile):
    """Write the slurm file"""
    # haoyuan:
    # This function set up the slurm file for one job
    fh = open(slurm_file_name,'w')
    fh.write("#!/bin/bash\n")
    fh.write("#SBATCH -p high2\n")
    fh.write("#SBATCH -n {:d}\n".format(ncpu))
    fh.write("#SBATCH --exclusive\n")
    fh.write("#SBATCH --mem=0\n")
    fh.write("#SBATCH --ntasks-per-node={:d}\n".format(tasks_per_node))
    #fh.write("#SBATCH -c 2\n")
    fh.write("#SBATCH --time=02:00:00\n")
    fh.write("#SBATCH --job-name={:s}\n".format(job_name))
    fh.write("#SBATCH --switches=1\n")
    fh.write("set -x\n")
    fh.write(environment_setup_commands + "\n")
    fh.write("module list\n")
    fh.write("ulimit -l unlimited\n")
    #fh.write("srun ./aspect {:s}\n".format(prmfile))
    #fh.write("mpirun -n {:d} --mca btl_openib_use_eager_rdma 1 --mca mpi_leave_pinned 1 --bind-to-core --report-bindings --mca btl_openib_allow_ib 1 ./aspect {:s}\n".format(ncpu,prmfile))
    #fh.write("mpirun -n {:d} --mca btl ^tcp --report-bindings ./aspect {:s}\n".format(ncpu,prmfile))
    # haoyuan: --bind-to core is always choosen as the default
    # one need to change the path to aspect to make this function
    # moreover, there could be a load library issure for openib
    # fh.write("mpirun -n {:d} --bind-to core --mca btl_openib_allow_ib 1 --mca btl openib,self,vader --report-bindings /home/lochy/software/aspect/build/aspect {:s}\n".format(ncpu,prmfile))
    fh.write("mpirun -n {:d} --bind-to core --report-bindings /home/lochy/software/aspect/build/aspect {:s}\n".format(ncpu,prmfile))
    #fh.write("mpirun -n {:d} ./aspect-impi {:s}\n".format(ncpu,prmfile))
    fh.close()
    
for core_count in core_counts:
    for resolution in refinement_levels:
        if( core_count >= minimum_core_count_for_refinement_level[resolution]
            and
            core_count <= maximum_core_count_for_refinement_level[resolution]):
            for setup in setups:
                jobname = "run_{:d}_{:d}_{:d}".format(core_count,resolution,setup)
                output_file = "tmp/{:s}/output_{:d}_{:d}_{:d}".format(cluster_label,core_count,resolution,setup)
                input_file = "tmp/{:s}/input_{:d}_{:d}_{:d}".format(cluster_label,core_count,resolution,setup)
                print(output_file)
                parameters = dict([])
                parameters['OUTPUT_DIRECTORY'] = output_file
                parameters['RESOLUTION'] = resolution
                
                # do string replacement on the base input file
                generate_input_file(base_input,input_file,parameters)
                slurm_file = input_file + ".slurm"
                # haoyuan: calls function to generate slurm file for one job
                generate_slurm_file(slurm_file,core_count,tasks_per_node,jobname,input_file)
                # haoyuan: submit the job
                os.system("sbatch " + slurm_file)
                

