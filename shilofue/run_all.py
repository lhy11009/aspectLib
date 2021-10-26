# -*- coding: utf-8 -*-
r"""Generates submission scripts (slurm-style) For the aspect performance benchmarks

haoyuan: this scripts test for different number of cores as well as different level of refinement
a tasks_per_node option is chosen, default is 32
make sure that you have the right commands in slurm before submitting them.

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python shilofue/run_all.py peloton-rome /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/
        -d master_TwoD

        -d: branch to test
  
  - hard in the code:

        bind-to options, in respective "generate_input_file*" functions
        
        core_counts, the number of cores you want to test
            
        openmpi, version of openmpi used
        
        module options, in respective "generate_input_file*" functions
    
        environment_setup_commands, set up mpi commands to use

        tasks_per_node option, in the 'main' function, governing the number of tasks per node
            
            for high2, use 64;
            for roma, use 128;

""" 

import numpy as np
import sys, argparse

# from subprocess import run
import os
import shutil

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
ASPECT_SOURCE_DIR = os.environ['ASPECT_SOURCE_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


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


def generate_slurm_file_peloton(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile,**kwargs):
    """Write the slurm file for peloton"""
    # modify this to contain the commands necessary to setup MPI environment
    #environment_setup_commands = "module load openmpi/3.1.3 intel-mkl"
    # haoyuan: I think we need to load 4.0.5 here
    environment_setup_commands="module unload openmpi/4.0.1\n\
            module load openmpi/4.1.0-mpi-io"
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
    # command to run
    branch=kwargs.get('branch', None)
    if branch==None:
        fh.write("mpirun -n {:d} --bind-to hwthread --report-bindings {:s}/build/aspect {:s}\n".format(ASPECT_SOURCE_DIR,ncpu,prmfile))
    else:
        fh.write("mpirun -n {:d} --bind-to hwthread --report-bindings {:s}/build_{:s}/aspect {:s}\n".format(ncpu, ASPECT_SOURCE_DIR,branch,prmfile))
    #fh.write("srun ./aspect {:s}\n".format(prmfile))
    #fh.write("mpirun -n {:d} --mca btl_openib_use_eager_rdma 1 --mca mpi_leave_pinned 1 --bind-to-core --report-bindings --mca btl_openib_allow_ib 1 ./aspect {:s}\n".format(ncpu,prmfile))
    #fh.write("mpirun -n {:d} --mca btl ^tcp --report-bindings ./aspect {:s}\n".format(ncpu,prmfile))
    
    # haoyuan: --bind-to core is always choosen as the default
    # one need to change the path to aspect to make this function
    # moreover, there could be a load library issure for openib
    # fh.write("mpirun -n {:d} --bind-to core --mca btl_openib_allow_ib 1 --mca btl openib,self,vader --report-bindings /home/lochy/software/aspect/build/aspect {:s}\n".format(ncpu,prmfile))
    # haoyuan: unload previous openmpi and reload a new one
    # fh.write("mpirun -n {:d} --bind-to hwthread --report-bindings /home/lochy/software/aspect/build/aspect {:s}\n".format(ncpu,prmfile))
    #fh.write("mpirun -n {:d} ./aspect-impi {:s}\n".format(ncpu,prmfile))
    fh.close()


def generate_slurm_file_peloton_rome(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile, **kwargs):
    """Write the slurm file for peloton roma
    
       Rome is the new partition for Magali's group 
    """
    # modify this to contain the commands necessary to setup MPI environment
    #environment_setup_commands = "module load openmpi/3.1.3 intel-mkl"
    # haoyuan: I think we need to load 4.0.5 here
    # haoyuan: unload previous openmpi and reload a new one
    environment_setup_commands="module unload openmpi/4.0.1\n\
            module load openmpi/4.1.0-mpi-io"
    # haoyuan:
    # This function set up the slurm file for one job
    fh = open(slurm_file_name,'w')
    fh.write("#!/bin/bash\n")
    fh.write("#SBATCH -p rome-256-512\n")
    fh.write("#SBATCH -n {:d}\n".format(ncpu))
    fh.write("#SBATCH --exclusive\n")
    # fh.write("#SBATCH --mem=0\n")
    fh.write("#SBATCH --ntasks-per-node={:d}\n".format(tasks_per_node))
    fh.write("#SBATCH --time=02:00:00\n")
    fh.write("#SBATCH --job-name={:s}\n".format(job_name))
    fh.write("#SBATCH --switches=1\n")
    fh.write("set -x\n")
    fh.write(environment_setup_commands + "\n")
    fh.write("module list\n")
    fh.write("ulimit -l unlimited\n")
    # haoyuan: --bind-to core is always choosen as the default
    # one need to change the path to aspect to make this function
    # moreover, there could be a load library issure for openib
    # command to run
    branch=kwargs.get('branch', None)
    if branch==None:
        fh.write("mpirun -n {:d} --bind-to socket --report-bindings {:s}/build/aspect {:s}\n".format(ncpu, ASPECT_SOURCE_DIR,prmfile))
    else:
        fh.write("mpirun -n {:d} --bind-to socket --report-bindings {:s}/build_{:s}/aspect {:s}\n".format(ncpu, ASPECT_SOURCE_DIR,branch,prmfile))
    # fh.write("mpirun -n {:d} --bind-to hwthread --report-bindings /home/lochy/software/aspect/build/aspect {:s}\n".format(ncpu,prmfile))
    fh.close()


def generate_slurm_file_stempede2(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile):
    """Write the slurm file"""
    # haoyuan:
    # This function set up the slurm file for one job
    fh = open(slurm_file_name,'w')
    fh.write("#!/bin/bash\n")
    fh.write("#SBATCH -p skx-normal\n")
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
    # haoyuan: unload previous openmpi and reload a new one
    # fh.write("module unload openmpi\nexport PATH=/home/rudolph/sw/openmpi-4.0.5/bin:$PATH\n")
    fh.write("mpirun -n {:d} --bind-to core --report-bindings /home1/06806/hyli/softwares/aspect {:s}\n".format(ncpu,prmfile))
    #fh.write("mpirun -n {:d} ./aspect-impi {:s}\n".format(ncpu,prmfile))
    fh.close()


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            server
        sys.arg[2](str):
            _path
    '''
    server = sys.argv[1] # server name
    _path = sys.argv[2] # path to the files
    # parse parameters
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str, default='', help='Some inputs')
    parser.add_argument('-o', '--outputs', type=str, default='.', help='Some outputs')
    parser.add_argument('-b', '--branch', type=str, default=None, help='The branch of code to test')
    _options = [] 
    try:
        _options = sys.argv[3: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)
    arg = parser.parse_args(_options)

    
    # Haoyuan: in this file, the field to change are marked with capital letters.
    # e.g.  set Output directory                       = OUTPUT_DIRECTORY
    # The 'base' input file that gets modified
    base_input = os.path.join(ASPECT_LAB_DIR, 'files', 'AffinityTest', 'spherical_shell_expensive_solver.prm')
    #cluster_label = "PI4CS_aspect-2.0-pre-40tasks"
    #cluster_label = "peloton-ii-64tasks-hwthread-openmpi-4.0.1"
    #cluster_label = "peloton-ii-32tasks-core-openmpi-4.0.1"
    
    
    
    # slurm parameterization
    # for peloton ii
    # core_counts = [1,2,4,8,16,32,64,128,256,512,768,1024]#,200,300,400]#,500,800,1000,1500]
    # for rome-256-512
    # 64 tasks per node
    core_counts = [1,2,4,8,16,32,64, 128,256, 512] # 768,1024]#,200,300,400]#,500,800,1000,1500]
    refinement_levels = [2,3,4,5]#,6]
    #                                          0   1   2   3       4     5    6
    minimum_core_count_for_refinement_level = [0,  0,   1,   1,   10, 100, 500]# for refinement levels 0-6
    maximum_core_count_for_refinement_level = [0,  0,1000,1000, 1000,2000,2000]
    
    setups = [1,]
    # for peloton ii
    # tasks_per_node = 32
    # for rome-256-512
    tasks_per_node = 128
    openmpi = "4.1.0"
    
    cluster_label = "%s-%stasks-socket-openmpi-%s" % (server, tasks_per_node, openmpi) # ?

    tmp_dir = os.path.join(_path, 'tmp')
    if not os.path.isdir(tmp_dir):
        # make a new directory every time
        # shutil.retree(tmp_dir)
        
        # make a new directory if old one doesn't exist
        os.mkdir(tmp_dir)

    cluster_dir = os.path.join(tmp_dir, cluster_label)
    if not os.path.isdir(cluster_dir):
        os.mkdir(cluster_dir)

    for core_count in core_counts:
        for resolution in refinement_levels:
            if( core_count >= minimum_core_count_for_refinement_level[resolution]
                and
                core_count <= maximum_core_count_for_refinement_level[resolution]):
                for setup in setups:
                    jobname = "run_{:d}_{:d}_{:d}".format(core_count,resolution,setup)
                    output_file = "{:s}/tmp/{:s}/output_{:d}_{:d}_{:d}".format(_path,cluster_label,core_count,resolution,setup)
                    input_file = "{:s}/tmp/{:s}/input_{:d}_{:d}_{:d}".format(_path, cluster_label,core_count,resolution,setup)
                    print(output_file)
                    parameters = dict([])
                    parameters['OUTPUT_DIRECTORY'] = output_file
                    parameters['RESOLUTION'] = resolution
                    
                    # do string replacement on the base input file
                    generate_input_file(base_input,input_file,parameters)
                    slurm_file = input_file + ".slurm"
                    # haoyuan: calls function to generate slurm file for one job
                    if server == "peloton-ii":
                        generate_slurm_file_peloton(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=arg.branch)
                        print("sbatch" + slurm_file)
                        os.system("sbatch " + slurm_file)
                    elif server == "stampede2":
                        generate_slurm_file(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=arg.branch)
                        print("sbatch" + slurm_file)
                        os.system("sbatch " + slurm_file)
                    if server == "peloton-rome":
                        generate_slurm_file_peloton_rome(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=arg.branch)
                        print("sbatch -A billen " + slurm_file)
                        os.system("sbatch -A billen " + slurm_file)
                    # haoyuan: submit the job
                    # os.system("sbatch " + slurm_file)


# run script
if __name__ == '__main__':
    main()
