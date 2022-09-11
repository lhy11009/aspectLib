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
import shilofue.ThDSubduction0.Cases as CasesThDSubduction
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
        python -m shilofue.AffinityTest create_tests -j `pwd`/affinity_test.json\n\
\n\
    analyze the results\n\
\n\
        python -m shilofue.AffinityTest analyze_results\
 -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_example\
 -c peloton-rome-128tasks-socket-openmpi-4.1.0 -o .\n\
        ")

    
def generate_input_file_1(base_file_name, output_file_name, output_dir, global_refinement, stokes_solver_type):
    """Read the 'base' input file from base_file_name, replace strings 
    using functions defined in ParsePrm.py, and write new output file to output_file_name"""
    with open(base_file_name,'r') as fh:
        idict = ParsePrm.ParseFromDealiiInput(fh)
    idict['Mesh refinement']["Initial global refinement"] = str(global_refinement)
    if "Stokes solver type" in idict['Solver parameters']['Stokes solver parameters'] or stokes_solver_type == "block GMG":
        # do this so that if the option is default, no need to create an entry
        idict['Solver parameters']['Stokes solver parameters']['Stokes solver type'] = stokes_solver_type
        if "Material averaging" not in idict["Material model"]:
            idict["Material model"]["Material averaging"] = "harmonic average"
    idict["Output directory"] = os.path.join("..", output_dir)
    output_dir = os.path.dirname(output_file_name)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)  # haoyuan: run module isn't imported correctly
    with open(output_file_name,'w') as ofh:
        ParsePrm.ParseToDealiiInput(ofh, idict)


class AFFINITY_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with AFFINITY
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        # todo_json
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Test directory", str, ["test directory"], ".", nick='test_dir')
        self.add_key("Base prm/json file (inputs)", str, ["base file"], "./test.prm", nick='base_file')
        self.add_key("Slurm file (inputs)", str, ["slurm file"], "./slurm.sh", nick='slurm_base_path')
        self.add_key("Server", str, ["server"], "peloton-high2", nick='server')
        self.add_key("Tasks per node", int, ["tasks per node"], 32, nick='tasks_per_node')
        self.add_key("Refinement level, note this is a summarized parameter of the refinement scheme assigned,\
it only takes effect if the input is positiveh",\
            list, ["refinement levels"], [1], nick="refinement_levels")
        self.add_key("Openmpi version", str, ["openmpi version"], "", nick='openmpi')
        self.add_key("Minimum cpus for each refinement level, repectively", list, ["minimum cores for refinement"], [], nick="min_cores_for_refinement")
        self.add_key("Maximum cpus for each refinement level, repectively", list, ["maximum cores for refinement"], [], nick="max_cores_for_refinement")
        self.add_key("project", str, ["project"], "", nick="project")
        self.add_key("branch", str, ["branch"], "master", nick="branch")
        self.add_key("List of nodes to test", list, ["node list"], [], nick="nodelist")
        self.add_key("End step", int, ["end step"], -1, nick="end_step")
        self.add_key("Stokes solver type", str, ["stokes solver type"], "block AMG", nick="stokes_type")
    
    def check(self):
        base_file = Utilities.var_subs(self.values[1])
        assert(os.path.isfile(base_file))
        assert(base_file.split('.')[-1] in ["json", "prm"])
        slurm_base_path = Utilities.var_subs(self.values[2])
        server = self.values[3]
        assert(server in ['peloton-rome', 'peloton-high2', 'stampede2'])
        refinement_levels = self.values[5]  # assert refinement levels are integers
        for refinement_level in refinement_levels:
            assert(type(refinement_level) == int) 
        min_cores_for_refinement = self.values[7]
        for core_count in min_cores_for_refinement:
            assert(type(core_count) == int)
        max_cores_for_refinement = self.values[8]
        for core_count in max_cores_for_refinement:
            assert(type(core_count) == int)
        stokes_type = self.values[13]
        assert(stokes_type in ["block AMG", "block GMG"])
    
    def to_init(self):
        '''
        interface to the init function of class AFFINITY
        '''
        test_dir = Utilities.var_subs(self.values[0])
        base_file = Utilities.var_subs(self.values[1])
        slurm_base_path = Utilities.var_subs(self.values[2])
        server = self.values[3]
        tasks_per_node = self.values[4]
        refinement_levels = self.values[5]
        end_step = self.values[12]
        stokes_type = self.values[13]
        return test_dir, base_file, slurm_base_path, server, tasks_per_node, refinement_levels, end_step, stokes_type
    
    def get_openmpi_version(self):
        openmpi = self.values[6]
        if openmpi != "":
            return openmpi
        else:
            return None
    
    def get_min_cores_for_refinement(self):
        min_cores_for_refinement = self.values[7]
        return min_cores_for_refinement
    
    def get_max_cores_for_refinement(self):
        max_cores_for_refinement = self.values[8]
        return max_cores_for_refinement
    
    def get_project(self):
        project = self.values[9]
        if project != "":
            return project
        else:
            return None
    
    def get_branch(self):
        branch = self.values[10]
        return branch
    
    def get_node_list(self):
        nodelist = self.values[11]
        return nodelist

    def get_test_dir(self):
        test_dir = Utilities.var_subs(self.values[0])
        return test_dir


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
    def __init__(self, test_dir, base_prm_path, slurm_base_path, server, tasks_per_node,\
                 refinement_levels, end_step, stokes_type, **kwargs):
        self.test_dir = test_dir
        self.base_prm_path = base_prm_path
        self.slurm_base_path = slurm_base_path
        openmpi = kwargs.get("openmpi", None)
        self.server = server
        self.tasks_per_node = tasks_per_node
        self.stokes_type = stokes_type
        self.cluster_label = "%s-%stasks-socket" % (self.server, tasks_per_node)
        if openmpi is not None:
            self.cluster_label += ("-openmpi-" + openmpi) # ?
        if self.stokes_type == "block GMG":
            self.cluster_label += "-bGMG"
        self.refinement_levels = refinement_levels
        self.end_step = end_step
        self.setups = [1, ] # unused
        self.core_counts = []
        for i in range(30):
            core_raw = 2.0**(i/2.0)
            if core_raw <= tasks_per_node:
                core = int(core_raw)
            else:
                core = int(core_raw//tasks_per_node * tasks_per_node)
            if not core in self.core_counts:
                self.core_counts.append(core)
        self.min_cores_for_refinement = kwargs.get("min_cores_for_refinement", [])
        if len(self.min_cores_for_refinement) > 0:
            assert(len(self.min_cores_for_refinement) == len(self.refinement_levels))
        else:
            self.min_cores_for_refinement = [0 for i in range(len(self.refinement_levels))]
        self.max_cores_for_refinement = kwargs.get("max_cores_for_refinement", [])
        if len(self.max_cores_for_refinement) > 0:
            assert(len(self.max_cores_for_refinement) == len(self.refinement_levels))
        else:
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
            elif self.project == "ThDSubduction":
                CASE = CasesThDSubduction.CASE
                CASE_OPT = CasesThDSubduction.CASE_OPT
            else:
                raise NotImplementedError("Options for project %s is not implemented." % self.project)
            # Call this in order to include all the related files
            CasesP.create_case_with_json(self.base_prm_path, CASE, CASE_OPT, reset_refinement_level=refinement,\
                fix_output_dir=os.path.dirname(input_dir), fix_case_name=os.path.basename(input_dir),\
                fix_case_output_dir="../%s" % os.path.basename(output_dir), end_step=self.end_step, reset_stokes_solver_type=self.stokes_type)
            print("case generated: %s" % prm_file)
        else:
            # There is only a prm file to take care of.
            generate_input_file_1(self.base_prm_path, prm_file, os.path.basename(output_dir), refinement, self.stokes_type)
            print("case generated: %s" % input_dir)
            # change the refinement level in the prm file
        # haoyuan: calls function to generate slurm file for one job
        SlurmOperator = ParsePrm.SLURM_OPERATOR(self.slurm_base_path)
        SlurmOperator.SetAffinity(np.ceil(core_count/self.tasks_per_node), core_count, 1)
        SlurmOperator.SetCommand(self.branch, os.path.basename(prm_file))
        SlurmOperator.SetName(jobname)
        SlurmOperator(slurm_file_output_path)
    
    def get_cluster_label(self):
        return self.cluster_label
    
    def get_test_inputs_dir(self):
        test_inputs_dir = "{:s}/tmp/{:s}".format(self.test_dir, self.cluster_label)
        return test_inputs_dir

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
        # print information of refinements and core counts
        for i in range(len(self.refinement_levels)):
            refinement = self.refinement_levels[i]
            print("refinement: ", refinement)
            is_first = True
            for core_count in self.core_counts:
                if( core_count >= self.min_cores_for_refinement[i]
                    and
                    core_count <= self.max_cores_for_refinement[i]):
                    if is_first:
                        print("core: ")
                        is_first = False
                    print(" %d" % core_count)
            print("\n")
        # generate test files
        for core_count in self.core_counts:
            for i in range(len(self.refinement_levels)):
                refinement = self.refinement_levels[i]
                if( core_count >= self.min_cores_for_refinement[i]
                    and
                    core_count <= self.max_cores_for_refinement[i]):
                    if is_first:
                        is_first = False
                    for setup in self.setups:
                        jobname = "run_{:d}_{:d}_{:d}".format(core_count,refinement,setup)
                        output_dir = "{:s}/tmp/{:s}/output_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        input_dir = "{:s}/tmp/{:s}/input_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        self.CreateCase(input_dir, output_dir, jobname, core_count, refinement)
                        # contents in the scripts: sbatch cases
                        bash_contents += "\ncd %s" % os.path.basename(input_dir)
                        bash_contents += "\nsbatch job.sh"
                        bash_contents += "\ncd .."
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
                prm_path = os.path.join(subdir, _dir, 'original.prm')
                target_path = os.path.join(results_dir, case_name + '.prm')
                print("copy %s to %s" % (prm_path, target_path))
                shutil.copy(prm_path, target_path)
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
        input_file = os.path.basename(output_file).split("_", maxsplit=1)[1] + ".prm"
        input_path = os.path.join(test_results_dir,  input_file)
        stokes_solver_type = ParsePrm.GetStokesSolverTypeFromPrm(input_path) 
        patterns = output_file.split('_')
        print("Output file found: %s" % output_file)
        # append data
        if stokes_solver_type == "block AMG":
            subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s" 
                            % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        elif stokes_solver_type == "block GMG":
            subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s gmg" 
                            % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        else:
            raise ValueError("stokes_solver_type must be block AMG or block GMG")
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


def create_tests_with_json(json_opt, AFFINITY, AFFINITY_OPT, **kwargs):
    '''
    A wrapper for the CASES class
    Inputs:
        json_opt(str, dict): path or dict a json file
        kwargs (dict):
            update (bool): update existing cases?
    Returns:
        case_dir: return case directory
    '''
    # todo_json
    Affinity_Opt = AFFINITY_OPT()
    if type(json_opt) == str:
        if not os.access(json_opt, os.R_OK):
            raise FileNotFoundError("%s doesn't exist" % json_opt)
        Affinity_Opt.read_json(json_opt)
    elif type(json_opt) == dict:
        Affinity_Opt.import_options(json_opt)
    else:
        raise TypeError("Type of json_opt must by str or dict")
    # check variables
    Affinity_Opt.check()
    test_dir = Affinity_Opt.get_test_dir()
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)
    Affinity = AFFINITY(*Affinity_Opt.to_init(), openmpi=Affinity_Opt.get_openmpi_version(),\
                        min_cores_for_refinement=Affinity_Opt.get_min_cores_for_refinement(),\
                        max_cores_for_refinement=Affinity_Opt.get_max_cores_for_refinement(),\
                        project=Affinity_Opt.get_project(), branch=Affinity_Opt.get_branch(),\
                        nodelist=Affinity_Opt.get_node_list())
    # remove older files
    inputs_dir = Affinity.get_test_inputs_dir()
    if os.path.isdir(inputs_dir):
        proceed = input("Target exists (%s), Delete older directory? (y/n)" % inputs_dir)  # check output dir exists
        if not proceed == 'y':
            exit(0)
        else:
            shutil.rmtree(inputs_dir)
    # generate test files
    Affinity()
    # save a copy of the json file
    json_output_path = os.path.join(Affinity_Opt.get_test_dir(),\
                                    "affinity_test_" + Affinity.get_cluster_label() + ".json")
    if type(json_opt) == str:
        shutil.copy2(json_opt, json_output_path)
    elif type(json_opt) == dict:
        with open(json_output_path, 'w') as fout:
            json.dump(fout, json_opt)
    print("Json file saved: %s" % json_output_path)

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
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='path to input file/directory')
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='path to an output file/directory')
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='path to a json file')
    parser.add_argument('-c', '--cluster', type=str,
                        default='',
                        help='name of the cluster')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if commend == 'create_tests':
        # todo_json
        create_tests_with_json(arg.json, AFFINITY, AFFINITY_OPT)

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
