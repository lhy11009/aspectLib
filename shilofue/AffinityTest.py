# -*- coding: utf-8 -*-
r"""Analyze affinity test results
The first part does the affinity test, while the second part does the analysis.
Notice that we might not have all the package installed on remote, so one should download the result and then run the analysis on there own laptop.

This outputs: 

  - a figure of results
""" 
import numpy as np
import sys, os, argparse
# on peloton,this is no such modules
import shutil
import numpy as np
import math
import shilofue.Cases as CasesP
import shilofue.TwoDSubduction0.Cases as CasesTwoDSubduction
try:
    import pathlib
    import subprocess
    from matplotlib import cm
    from matplotlib import pyplot as plt
except ImportError:
    pass

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
ASPECT_SOURCE_DIR = os.environ['ASPECT_SOURCE_DIR']


def Usage():
    print("\
Run affinity test and analyze results\n\
The first part does the affinity test, while the second part does the analysis.\n\
Notice that we might not have all the package installed on remote, so one should download the result and then run the analysis on there own laptop.\n\
\n\
Examples of usage: \n\
\n\
  - Do affinity test: \n\
\n\
    (run this part on server)\n\
\n\
    example:\n\
\n\
        (substitute -i and -o options with your own prm file and output directory)\n\
        python shilofue/AffinityTest.py run_tests -s peloton-rome -t 128\
 -i /home/lochy/ASPECT_PROJECT/aspectLib/files/AffinityTest/spherical_shell_expensive_solver.prm\
 -o $TwoDSubduction_DIR/rene_affinity_test\
 -nl topaz-0 topaz-2 topaz-3\
 -m 385\n\
\n\
  - Analyze results: \n\
    (run this part on a laptop)\n\
\n\
    example:\n\
\n\
        python shilofue/AffinityTest.py analyze_results\
 -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_example\
 -c peloton-rome-128tasks-socket-openmpi-4.1.0 -o .\n\
        ")

    
def my_assert(_condition, _errortype, _message):
    '''
    an assert function for runtime use
    Inputs:
        _condition(True or False):
            the condition to assert
        _errortype(some error):
            the type of error to raise
        _message(str):
            the message to raise
    '''
    if _condition is False:
        raise _errortype(_message)


def generate_input_file(base_file_name,output_file_name,dictionary):
    """Read the 'base' input file from base_file_name, replace strings 
    using dictionary, and write new output file to output_file_name"""
    fh = open(base_file_name,'r')
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)  # haoyuan: run module isn't imported correctly
    ofh = open(output_file_name,'w')
    for line in fh:        
        for key in dictionary:
            if key in line:                
                line = line.replace(key,str(dictionary[key]))
        ofh.write(line)
    fh.close()
    ofh.close()


def generate_slurm_file_peloton(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile,**kwargs):
    """
    Write the slurm file for peloton.
    Inputs:
        slurm_file_name(str): full path for the slurm file
        ncpu (int): numbers of cpu to use
        tasks_per_node (int): tasks per node
        job_name (str)
        prmfile: full path for the prm file
        kwargs(dict):
            branch (str): the branch of aspect to run with
    """
    # modify this to contain the commands necessary to setup MPI environment
    # environment_setup_commands = "module load openmpi/3.1.3 intel-mkl"
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


def generate_slurm_file_peloton_roma(slurm_file_name,ncpu,tasks_per_node,job_name,prmfile, **kwargs):
    """Write the slurm file for peloton roma
        Inputs:
            slurm_file_name(str): full path for the slurm file
            ncpu (int): numbers of cpu to use
            tasks_per_node (int): tasks per node
            job_name (str)
            prmfile: full path for the prm file
            kwargs(dict):
                branch (str): the branch of aspect to run with
                nodelist: (list of str) - list of node to run on
    """
    # modify this to contain the commands necessary to setup MPI environment
    #environment_setup_commands = "module load openmpi/3.1.3 intel-mkl"
    # haoyuan: I think we need to load 4.0.5 here
    # haoyuan: unload previous openmpi and reload a new one
    environment_setup_commands="module unload openmpi/4.0.1\n\
            module load openmpi/4.1.0-mpi-io"
    nodelist= kwargs.get('nodelist', [])
    assert(type(ncpu) == int)
    assert(type(tasks_per_node) == int)
    assert(type(nodelist) == list)
    optional = ""
    # append nodelist
    if len(nodelist) > 0:
        nnode = int(math.ceil(1.0*ncpu/tasks_per_node))  # number of nodes required by size
        optional += "#SBATCH --nodelist="
        for i in range(min(len(nodelist), nnode)):
            node = nodelist[i]
            if i == 0:
                optional += node
            else:
                optional += ",%s" % node
        optional += "\n"
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
    if optional != "":
        fh.write(optional) # extra configuration
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
    """Write the slurm file for stempede2
        Inputs:
            slurm_file_name(str): full path for the slurm file
            ncpu (int): numbers of cpu to use
            tasks_per_node (int): tasks per node
            job_name (str)
            prmfile: full path for the prm file
            kwargs(dict):
                branch (str): the branch of aspect to run with
    """
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


def do_tests(server, _path, tasks_per_node, base_input_path, **kwargs):
    '''
    perform the affinity test
    Inputs:
        server(str): options for server
        _path(str): directory to hold test files, pull path
        tasks_per_node (int): number of tasks to run on each node.
        base_input_path: The 'base' input file that gets modified
        kwargs(dict):
            branch (str): branch to use, default is master
            nodelist: (list of str) - list of node to run on
            debug (int) : debug mode (1), normal (0)
            max_core_count (int): maximum number for core count
    '''
    # Haoyuan: in this file, the field to change are marked with capital letters.
    # e.g.  set Output directory                       = OUTPUT_DIRECTORY
    nodelist= kwargs.get('nodelist', [])  # list of nodes
    assert(type(nodelist) == list)
    if len(nodelist) > 0:
        assert(server in ["peloton-rome"]) # this only works with these servers
    max_core_count = kwargs.get('max_core_count', 10000)
    debug = kwargs.get('debug', 0) # debug mode
    project = kwargs.get('project', None)  # use a project instead of the default prm
    json_file = kwargs.get('json_file', None)
    # slurm parameterization
    # for peloton ii
    # core_counts = [1,2,4,8,16,32,64,128,256,512,768,1024]#,200,300,400]#,500,800,1000,1500]
    # for rome-256-512
    # 64 tasks per node
    setups = [1, ]
    all_available_core_counts = [1,2,4,8,16,32,64, 128,256, 512, 768,1024]#,200,300,400]#,500,800,1000,1500]
    if debug:
        core_counts = [1,2,4] # debug mode, see all the settings work
    else:
        core_counts = []
        for core_count in all_available_core_counts:
            if core_count < max_core_count:
                core_counts.append(core_count)
    refinement_levels = [2,3,4,5]#,6]
    # refinement_levels = [2]
    #                                          0   1   2   3       4     5    6
    minimum_core_count_for_refinement_level = [0,  0,   1,   1,   10, 100, 500]# for refinement levels 0-6
    maximum_core_count_for_refinement_level = [0,  0,1000,1000, 1000,2000,2000] 
    openmpi = "4.1.0"  # change with updates
    cluster_label = "%s-%stasks-socket-openmpi-%s" % (server, tasks_per_node, openmpi) # ?
    if len(nodelist) > 0:
        for i in range(len(nodelist)):
            node = nodelist[i]
            if i == 0:
                cluster_label += "-%s" % node
            else:
                cluster_label += "_%s" % node
    # create a directory to hold test results
    tmp_dir = os.path.join(_path, 'tmp')
    if not os.path.isdir(tmp_dir):
        # make a new directory every time
        # shutil.retree(tmp_dir)
        # make a new directory if old one doesn't exist
        os.mkdir(tmp_dir)
    # make a subdirectory with the name of the cluster
    cluster_dir = os.path.join(tmp_dir, cluster_label)
    if not os.path.isdir(cluster_dir):
        os.mkdir(cluster_dir)
    # generate test files
    branch = kwargs.get('branch', 'master')
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
                    # todo_affinity
                    if project is not None:
                        assert(os.path.isfile(json_file))
                        if project == "TwoDSubduction":
                            CASE = CasesTwoDSubduction.CASE
                            CASE_OPT = CasesTwoDSubduction.CASE_OPT
                        else:
                            NotImplementedError("Options for project %s is not implemented." % project)
                        CasesP.create_case_with_json(json_file, CASE, CASE_OPT, reset_refinement_level=parameters["RESOLUTION"],\
                            fix_output_dir=parameters["OUTPUT_DIRECTORY"])
                    else:
                        generate_input_file(base_input_path,input_file,parameters)
                    slurm_file = input_file + ".slurm"
                    # haoyuan: calls function to generate slurm file for one job
                    if server == "peloton-ii":
                        generate_slurm_file_peloton(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=branch)
                        command_line="sbatch " + slurm_file
                    elif server == "stampede2":
                        generate_slurm_file(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=branch)
                        command_line="sbatch " + slurm_file
                    if server == "peloton-rome":
                        generate_slurm_file_peloton_roma(slurm_file,core_count,tasks_per_node,jobname,input_file, branch=branch, nodelist=nodelist)
                        command_line="sbatch " + "-A billen " + slurm_file
                    print(command_line)
                    os.system(command_line)

def organize_result(test_root_dir, cluster):
    '''
    organize result
    '''
    # check directory existance
    _dir = os.path.join(test_root_dir, 'results')
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    _dir = os.path.join(test_root_dir, 'results', 'spherical_shell_expensive_solver')
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    _dir = os.path.join(test_root_dir, 'results', 'spherical_shell_expensive_solver', cluster)
    if not os.path.isdir(_dir):
        os.mkdir(_dir)
    results_dir = _dir
    results_tmp_dir = os.path.join(test_root_dir, 'tmp', cluster)
    assert(os.path.isdir(results_tmp_dir))
    for subdir, dirs, _ in os.walk(results_tmp_dir):
        for _dir in dirs:
            if _dir.startswith('output'):
                case_name = _dir.split('_', 1)[1]
                log_path = os.path.join(subdir, _dir, 'log.txt')
                target_path = os.path.join(results_dir, 'output_' + case_name)
                print("copy %s to %s" % (log_path, target_path))
                shutil.copy(log_path, target_path)
            else:
                continue

    return results_dir


def analyze_affinity_test_results(test_results_dir, output_dir):
    '''
    analyze affinity test results
    '''
    total_wall_clock = []
    assemble_stokes_system = []
    solve_stokes_system = []
    cores = []
    resolutions = []
    setups = []
    # go into sub dirs
    temp_file = os.path.join(ASPECT_LAB_DIR, 'temp')  # file to save partial results
    path_obj = pathlib.Path(test_results_dir).rglob("output*")
    i = 0
    for _path in path_obj:
        i += 1
        output_file = str(_path)
        output_path = os.path.join(test_results_dir, output_file)
        patterns = output_file.split('_')
        print("Output file found: %s" % output_path)
        # append data
        subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s" 
                       % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        try:
            data = np.genfromtxt(temp_file)
            total_wall_clock.append(data[0, -1])
            assemble_stokes_system.append(data[1, -1])
            solve_stokes_system.append(data[2, -1])
        except Exception:
            pass
        else:
            setups.append(int(patterns[-1]))
            resolutions.append(int(patterns[-2]))
            cores.append(int(patterns[-3]))
    my_assert(i > 0, AssertionError, "There is no output* file in the folder %s" % test_results_dir)

    setups = np.array(setups)
    resolutions = np.array(resolutions)
    cores = np.array(cores)
    total_wall_clock = np.array(total_wall_clock)
    assemble_stokes_system = np.array(assemble_stokes_system)
    solve_stokes_system = np.array(solve_stokes_system)
    # rearrange and sort data
    print("Affinity Test Results:")
    print('setups:', setups)
    print('resolutions', resolutions)
    print("cores:", cores)
    print("total wall clocks:", total_wall_clock)
    print("assemble_stokes_system:", assemble_stokes_system)
    print("solve_stokes_system:", solve_stokes_system)
    
    # plot via matplotlib
    resolution_options = []
    for resolution in resolutions:
        if resolution not in resolution_options:
            resolution_options.append(resolution)
    resolution_options = np.array(resolution_options)
    resolution_options = np.sort(resolution_options)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        plot_indexes = (resolutions == resolution) 
        xs = cores[plot_indexes]
        sequential_indexes = np.argsort(xs)
        ys1 = total_wall_clock[plot_indexes]
        # labels
        _label0 = 'Total Wall Clock(resolution=%d)' % resolution
        _label1 = 'Assemble Stokes System(resolution=%d)' % resolution
        _label2 = 'Solve Stokes System(resolution=%d)' % resolution
        # plot
        axs[0].loglog(xs[sequential_indexes], ys1[sequential_indexes], ".-", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label0)
        ys2 = assemble_stokes_system[plot_indexes]
        axs[0].loglog(xs[sequential_indexes], ys2[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label1)
        ys3 = ys2 / ys1
        axs[1].semilogx(xs[sequential_indexes], ys3[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label1)
        ys4 = solve_stokes_system[plot_indexes]
        axs[0].loglog(xs[sequential_indexes], ys4[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label2)
        ys5 = ys4 / ys1
        axs[1].semilogx(xs[sequential_indexes], ys5[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label=_label2)
    axs[0].set_xlabel('Cores')
    axs[0].set_ylabel('Time [s]')
    axs[0].grid()
    axs[0].set_title('Wall Clock')
    axs[0].legend(fontsize='x-small')
    axs[1].set_xlabel('Cores')
    axs[1].set_ylabel('Percentage')
    axs[1].grid()
    axs[1].set_title('Percentage of Each Part')
    # title and save path 
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    plt.savefig(filepath)

    pass

def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--inputs', type=str, default='.',\
            help='Path of input.\n\
            In case of running the test, this is the path to the prm file\n\
            In case of analyze the results, this is the path of the directory for testing results')
    parser.add_argument('-o', '--outputs', type=str, default='.',\
            help='Path of output.\n\
            In case of running the test, this is the place to save results\n\
            In case of analyze the results, this is the place to save images')
    parser.add_argument('-c', '--cluster', type=str, default='peloton-ii-32tasks-core-openmpi-4.0.1',\
            help='name of the cluster.\n\
            This information is only needed when doing analysis.\n\
            At that time, you need to look into the data directory for the testing results.\n\
            This variable is the name of the subdirectory in there')
    parser.add_argument('-s', '--server', type=str, default='peloton-rome',\
            help='server and partition to run tests on.\n\
            one in (peloton-ii, peloton-rome, stampede2)')
    parser.add_argument('-b', '--branch', type=str, default=None, help='The branch of aspect to test')
    parser.add_argument('-t', '--tasks_per_node', type=int, default=32, help='Number of task to run on each node')
    parser.add_argument('-nl', '--nodelist', nargs="+", default=[], help='list of nodes to run on, separate with blank in inputs')
    parser.add_argument('-d', '--debug', type=int, default=0,\
            help='Run in debug mode by only deploying\
            a few cases with low number of nodes')
    parser.add_argument('-m', '--max_cpu', type=int, default=1e10,\
            help="Maximum number of cpu assigned.\
            This is not the real maximum number of cpus to be run on,\
            but meant to manually put a limit on that.\
            This doesn't work with the debug mode.")
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if commend == 'run_tests':
        # base_input_path = os.path.join(ASPECT_LAB_DIR, 'files', 'AffinityTest', 'spherical_shell_expensive_solver.prm')
        assert(os.path.isdir(arg.outputs))  # check output dir exists
        do_tests(arg.server, arg.outputs, arg.tasks_per_node, arg.inputs, branch=arg.branch, nodelist=arg.nodelist, debug=arg.debug, max_core_count=arg.max_cpu)  # create and submit jobs
    elif commend == 'analyze_results':
        # example:
        # python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        # -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_20211025 -c peloton-rome-128tasks-socket-openmpi-4.1.0
        # todo
        assert(os.path.isdir(arg.outputs))  # check output dir exists
        test_results_dir = organize_result(arg.inputs, arg.cluster)
        analyze_affinity_test_results(test_results_dir, arg.outputs)
    elif commend in ['help', '-h']:
        Usage()
        arg = parser.parse_args(['-h'])
    else:
        raise ValueError("No such command (%s), see Usage (run with help or -h)" % commend)

# run script
if __name__ == '__main__':
    main()
