# -*- coding: utf-8 -*-
r"""Analyze affinity test results
The first part does the affinity test, while the second part does the analysis.
Notice that we might not have all the package installed on remote, so one should download the result and then run the analysis on there own laptop.

This outputs: 

  - a figure of results
""" 
from sre_constants import SRE_FLAG_UNICODE
import numpy as np
import sys, os, argparse
# on peloton,this is no such modules
import shutil
import numpy as np
import math
import shilofue.Cases as CasesP
import shilofue.TwoDSubduction0.Cases as CasesTwoDSubduction
import shilofue.ParsePrm as ParsePrm
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

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


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
        (-t option gives the number of tasks per node, check this with the set ups of server)\n\
        (-m is a manual upper limit of the cpus to use)\n\
        (-rf is a list of refinment level to test)\n\
        (minc and maxc are lists of the minimum cores and maximum cores to run on, for each level of refinement, respectively.)\n\
        python -m shilofue.AffinityTest run_tests -s peloton-rome -t 128\
 -i /home/lochy/ASPECT_PROJECT/aspectLib/files/AffinityTest/spherical_shell_expensive_solver.prm\
 -o $TwoDSubduction_DIR/rene_affinity_test\
 -nl topaz-0 topaz-2 topaz-3\
 -m 385\
 -rf 2 3\
 -minc 1 1 -maxc 16 16\n\
\n\
        (for running affinity tests for a research project: -i is the json for the project cases instead, -p is the name\
of the project.\n\
        Note on doing this for a new project: search for \"project\" in this script and apply changes)\n\
        python -m shilofue.AffinityTest run_tests -i /home/lochy/ASPECT_PROJECT/aspectLib/files/TwoDSubduction/220708/case.json\
 -o /home/lochy/ASPECT_PROJECT/aspectLib/results -p TwoDSubduction\n\
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

    
def generate_input_file_1(base_file_name, output_file_name, output_path, global_refinement):
    """Read the 'base' input file from base_file_name, replace strings 
    using functions defined in ParsePrm.py, and write new output file to output_file_name"""
    with open(base_file_name,'r') as fh:
        idict = ParsePrm.ParseFromDealiiInput(fh)
    idict['Mesh refinement']["Initial global refinement"] = str(global_refinement)
    idict["Output directory"] = output_path
    output_dir = os.path.dirname(output_file_name)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)  # haoyuan: run module isn't imported correctly
    with open(output_file_name,'w') as ofh:
        ParsePrm.ParseToDealiiInput(ofh, idict)


class AFFINITY():
    '''
    class for running affinity tests
    Attributes:
        server(str): options for server
        test_dir(str): directory to hold test files, pull path
        tasks_per_node (int): number of tasks to run on each node.
        base_prm_path: The 'base' input file that gets modified, if project is None,
            this option points to a json file to import the settings
        branch (str): branch to use, default is master
        nodelist: (list of str) - list of node to run on
        debug (int) : debug mode (1), normal (0)
        max_core_count (int): maximum number for core count
    '''
    def __init__(self, test_dir, base_prm_path, slurm_base_path, server, tasks_per_node, refinement_levels, **kwargs):
        self.test_dir = test_dir
        self.base_prm_path = base_prm_path
        self.slurm_base_path = slurm_base_path
        self.core_counts = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 768, 1024]
        openmpi = kwargs.get("openmpi", None)
        self.server = server
        self.tasks_per_node = tasks_per_node
        self.cluster_label = "%s-%stasks-socket" % (self.server, tasks_per_node)
        self.cluster_label += ("-openmpi-" + openmpi) # ?
        self.refinement_levels = refinement_levels
        self.setups = [1, ] # unused
        try:
            self.min_cores_for_refinement = kwargs["min_cores_for_refinement"]
            assert(len(self.min_cores_for_refinement) == len(self.refinement_levels))
        except KeyError:
            self.min_cores_for_refinement = [0 for i in range(len(self.refinement_levels))]
        try:
            self.max_cores_for_refinement = kwargs["max_cores_for_refinement"]
            assert(len(self.max_cores_for_refinement) == len(self.refinement_levels))
        except KeyError:
            self.max_cores_for_refinement = [1e31 for i in range(len(self.refinement_levels))]
        self.project = kwargs.get("project", None)
        self.branch = kwargs.get("branch", "master")
        self.nodelist= kwargs.get('nodelist', [])  # list of nodes
        pass

    def CreateCase(self, input_dir, output_dir, jobname, core_count, refinement):
        '''
        create one case
        '''               
        slurm_file_output_path = os.path.join(input_dir, "job.sh")
        prm_file = os.path.join(input_dir, "case.prm")
        if not os.path.isdir(os.path.dirname(input_dir)):
            os.mkdir(os.path.dirname(input_dir))
        
        # do string replacement on the base input file
        if self.project is not None:
            assert(os.path.isfile(self.base_prm_path))
            if self.project == "TwoDSubduction":
                CASE = CasesTwoDSubduction.CASE
                CASE_OPT = CasesTwoDSubduction.CASE_OPT
            else:
                NotImplementedError("Options for project %s is not implemented." % self.project)
            CasesP.create_case_with_json(self.base_prm_path, CASE, CASE_OPT, reset_refinement_level=refinement,\
                fix_output_dir=os.path.dirname(input_dir), fix_case_name=os.path.basename(input_dir),\
                fix_case_output_dir="../%s" % os.path.basename(output_dir))
            print("case generated: %s" % prm_file)
        else:
            generate_input_file_1(self.base_prm_path, prm_file, os.path.basename(output_dir), refinement)
            print("case generated: %s" % input_dir)
            # change the refinement level in the prm file
        # haoyuan: calls function to generate slurm file for one job
        SlurmOperator = ParsePrm.SLURM_OPERATOR(self.slurm_base_path)
        SlurmOperator.SetAffinity(np.ceil(core_count/self.tasks_per_node), core_count, 1)
        SlurmOperator.SetCommand(self.branch, os.path.basename(prm_file))
        SlurmOperator.SetName(jobname)
        SlurmOperator(slurm_file_output_path)

    def __call__(self):
        # create a directory to hold test results
        tmp_dir = os.path.join(self.test_dir, 'tmp')
        if not os.path.isdir(tmp_dir):
            # make a new directory every time
            # shutil.retree(tmp_dir)
            # make a new directory if old one doesn't exist
            os.mkdir(tmp_dir)
        # make a subdirectory with the name of the cluster
        cluster_dir = os.path.join(tmp_dir, self.cluster_label)
        if not os.path.isdir(cluster_dir):
            os.mkdir(cluster_dir)
        # scripts outputs
        bash_output_path = os.path.join(cluster_dir, 'run_tests.sh')
        bash_contents = "#!/bin/bash" # scripts for running all the tests
        # generate test files
        print(self.core_counts)  # debug
        print(self.refinement_levels)
        for core_count in self.core_counts:
            for i in range(len(self.refinement_levels)):
                refinement = self.refinement_levels[i]
                if( core_count >= self.min_cores_for_refinement[i]
                    and
                    core_count <= self.max_cores_for_refinement[i]):
                    for setup in self.setups:
                        jobname = "run_{:d}_{:d}_{:d}".format(core_count,refinement,setup)
                        output_dir = "{:s}/tmp/{:s}/output_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        input_dir = "{:s}/tmp/{:s}/input_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        self.CreateCase(input_dir, output_dir, jobname, core_count, refinement)
                        # contents in the scripts: sbatch cases
                        bash_contents += "\ncd %s" % os.path.basename(input_dir)
                        bash_contents += "\nsbatch job.sh"
                        if self.project is not None:
                            bash_contents += "\ncd .."
        with open(bash_output_path, 'w') as fout:
            fout.write(bash_contents)


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
    Utilities.my_assert(i > 0, AssertionError, "There is no output* file in the folder %s" % test_results_dir)

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
    parser.add_argument('-rf', '--refinement_levels', nargs="+", default=[], help='list of refinement level to run on, separate with blank in inputs')
    parser.add_argument('-minc', '--min_cores_for_refinement', nargs="+", default=[], help='list of minimum cores to run on, separate with blank in inputs')
    parser.add_argument('-maxc', '--max_cores_for_refinement', nargs="+", default=[], help='list of maximum cores to run on, separate with blank in inputs')
    parser.add_argument('-d', '--debug', type=int, default=0,\
            help='Run in debug mode by only deploying\
            a few cases with low number of nodes')
    parser.add_argument('-m', '--max_cpu', type=int, default=1e10,\
            help="Maximum number of cpu assigned.\
            This is not the real maximum number of cpus to be run on,\
            but meant to manually put a limit on that.\
            This doesn't work with the debug mode.")
    parser.add_argument('-p', '--project', type=str, default=None, help='A specific project to test, for default tests, leave this blank')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if commend == 'run_tests':
        assert(os.path.isdir(arg.outputs))  # check output dir exists
        # do_tests(arg.server, arg.outputs, arg.tasks_per_node, arg.inputs, branch=arg.branch, nodelist=arg.nodelist,\
        #    debug=arg.debug, max_core_count=arg.max_cpu, project=arg.project)  # create and submit jobs
        openmpi =  "4.1.0"
        branch = "master_TwoD"
        slurm_base_path = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "220810", "job_p-billen.sh")
        refinement_levels = [int(rf) for rf in arg.refinement_levels]
        min_cores_for_refinement = [int(mc) for mc in arg.min_cores_for_refinement]
        max_cores_for_refinement = [int(mc) for mc in arg.max_cores_for_refinement]
        Affinity = AFFINITY(arg.outputs, arg.inputs, slurm_base_path, arg.server, arg.tasks_per_node,\
            refinement_levels, openmpi=openmpi, branch=branch, min_cores_for_refinement=min_cores_for_refinement,\
            max_cores_for_refinement=max_cores_for_refinement, project=arg.project)
        Affinity()

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
