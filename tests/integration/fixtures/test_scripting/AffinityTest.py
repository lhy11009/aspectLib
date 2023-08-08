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
from shilofue.PlotStatistics import STATISTICS_PLOT
from matplotlib import gridspec
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
        Lib_AffinityTest create_tests -j `pwd`/affinity_test.json\n\
\n\
    analyze the results\n\
\n\
        Lib_AffinityTest analyze_results\
 -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_example\
 -c peloton-rome-128tasks-socket-openmpi-4.1.0\n\
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
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Test directory", str, ["test directory"], ".", nick='test_dir')
        self.add_key("Base prm/json file (inputs)", list, ["base files"], [], nick='base_files')
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
        self.add_key("build directory", str, ["build directory"], "", nick="build_directory")
        self.add_key("List of nodes to test", list, ["node list"], [], nick="nodelist")
        self.add_key("End step", int, ["end step"], -1, nick="end_step")
        self.add_key("Stokes solver type", str, ["stokes solver type"], "block AMG", nick="stokes_type")
        self.add_key("Flag", str, ["flag"], "", nick="flag")
        self.add_key("Use mpirun", int, ["use mpirun"], 0, nick="use_mpirun")
        self.add_key("Bind to option", str, ["bind to"], "", nick="bind_to")
    
    def check(self):
        base_files = self.values[1]
        assert(len(base_files) > 0)
        for base_file in base_files:
            # assert file exist
            # assert file format is either json or prm
            base_file = Utilities.var_subs(base_file)
            Utilities.my_assert(FileExistsError, os.path.isfile(base_file),\
                                "%s is not a file" % base_file)
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
        bind_to = self.values[16]
        assert(bind_to in ["", "socket", "core"])
    
    def to_init(self):
        '''
        interface to the init function of class AFFINITY
        '''
        test_dir = Utilities.var_subs(self.values[0])
        # substitue base file path with absolute path
        base_files = []
        base_files_raw = self.values[1]
        for base_file in base_files_raw:
            base_file = Utilities.var_subs(base_file)
            base_files.append(base_file)
        slurm_base_path = Utilities.var_subs(self.values[2])
        server = self.values[3]
        tasks_per_node = self.values[4]
        refinement_levels = self.values[5]
        end_step = self.values[12]
        stokes_type = self.values[13]
        flag = self.values[14]
        use_mpirun = self.values[15]
        bind_to = self.values[16]
        return test_dir, base_files, slurm_base_path, server, tasks_per_node, refinement_levels, end_step,\
        stokes_type, flag, use_mpirun, bind_to
    
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
    
    def get_build_directory(self):
        build_directory = self.values[10]
        return build_directory
    
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
        base_files: The 'base' input file that gets modified, if project is None,
            this option points to a json file to import the settings
        build_directory (str): build_directory to use, default is master
        nodelist: (list of str) - list of node to run on
        debug (int) : debug mode (1), normal (0)
        max_core_count (int): maximum number for core count
    '''
    def __init__(self, test_dir, base_files, slurm_base_path, server, tasks_per_node,\
                 refinement_levels, end_step, stokes_type, flag, use_mpirun, bind_to, **kwargs):
        self.test_dir = test_dir
        self.base_files = base_files
        self.slurm_base_path = slurm_base_path
        openmpi = kwargs.get("openmpi", None)
        self.server = server
        self.tasks_per_node = tasks_per_node
        self.stokes_type = stokes_type
        self.cluster_label = "%s-%stasks-socket" % (self.server, tasks_per_node)
        if openmpi is not None:
            self.cluster_label += ("-openmpi-" + openmpi) # add version of openmpi
        if self.stokes_type == "block GMG":
            self.cluster_label += "-bGMG"  # add the type of Stokes solver
        self.refinement_levels = refinement_levels
        self.end_step = end_step
        # a list from the number of base files
        self.setups = [i+1 for i in range(len(base_files))]
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
        self.build_directory = kwargs.get("build_directory", "")
        self.nodelist= kwargs.get('nodelist', [])  # list of nodes
        self.use_mpirun = use_mpirun
        if self.use_mpirun:
            self.cluster_label += "-mpirun"
        if bind_to == "":
            self.bind_to = None
        else:
            self.bind_to = bind_to
            self.cluster_label += "-bind_to_%s" % self.bind_to
        if flag != "":
            self.cluster_label += "-%s" % flag  # add an assigned flag
        pass

    def CreateCase(self, setup, input_dir, output_dir, jobname, core_count, refinement):
        '''
        create one case
        Inputs:
            setup: the index of the input file to use
        '''               
        slurm_file_output_path = os.path.join(input_dir, "job.sh")
        prm_file = os.path.join(input_dir, "case.prm")
        if not os.path.isdir(os.path.dirname(input_dir)):
            os.mkdir(os.path.dirname(input_dir))
        
        # do string replacement on the base input file
        base_file  = self.base_files[setup-1]
        Utilities.my_assert(FileExistsError, os.path.isfile(base_file),\
                            "%s is not a file." % base_file)
        if self.project is not None:
            if self.project == "TwoDSubduction":
                CASE = CasesTwoDSubduction.CASE
                CASE_OPT = CasesTwoDSubduction.CASE_OPT
            elif self.project == "ThDSubduction":
                CASE = CasesThDSubduction.CASE
                CASE_OPT = CasesThDSubduction.CASE_OPT
            else:
                raise NotImplementedError("Options for project %s is not implemented." % self.project)
            # Call this in order to include all the related files
            CasesP.create_case_with_json(base_file, CASE, CASE_OPT, reset_refinement_level=refinement,\
                fix_output_dir=os.path.dirname(input_dir), fix_case_name=os.path.basename(input_dir),\
                fix_case_output_dir="../%s" % os.path.basename(output_dir), end_step=self.end_step, reset_stokes_solver_type=self.stokes_type)
            print("case generated: %s" % prm_file)
        else:
            # There are only prm files to take care of (base files are prms)
            generate_input_file_1(base_file, prm_file, os.path.basename(output_dir), refinement, self.stokes_type)
            print("case generated: %s" % input_dir)
        # haoyuan: calls function to generate slurm file for one job
        SlurmOperator = ParsePrm.SLURM_OPERATOR(self.slurm_base_path)
        SlurmOperator.SetAffinity(np.ceil(core_count/self.tasks_per_node), core_count, 1, use_mpirun=self.use_mpirun, bind_to=self.bind_to)
        SlurmOperator.SetCommand(self.build_directory, os.path.basename(prm_file))
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
        for setup in self.setups:
            for i in range(len(self.refinement_levels)):
                refinement = self.refinement_levels[i]
                for core_count in self.core_counts:
                    if( core_count >= self.min_cores_for_refinement[i]
                        and
                        core_count <= self.max_cores_for_refinement[i]):
                        if is_first:
                            is_first = False
                        jobname = "run_{:d}_{:d}_{:d}".format(core_count,refinement,setup)
                        output_dir = "{:s}/tmp/{:s}/output_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        input_dir = "{:s}/tmp/{:s}/input_{:d}_{:d}_{:d}".format(self.test_dir, self.cluster_label,core_count,refinement,setup)
                        self.CreateCase(setup, input_dir, output_dir, jobname, core_count, refinement)
                        # contents in the scripts: sbatch cases
                        bash_contents += "\ncd %s" % os.path.basename(input_dir)
                        bash_contents += "\nsbatch job.sh"
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
    # copy & paste results
    # 1. log outputs
    # 2. prm outputs
    # 3. statistics outputs
    assert(os.path.isdir(results_tmp_dir))
    for _dir in os.listdir(results_tmp_dir):
        if _dir.startswith('output'):
            case_name = _dir.split('_', 1)[1]
            log_path = os.path.join(results_tmp_dir, _dir, 'log.txt')
            assert(os.path.isfile(log_path))
            target_path = os.path.join(results_dir, 'output_' + case_name)
            print("copy %s to %s" % (log_path, target_path))
            shutil.copy(log_path, target_path)
            prm_path = os.path.join(results_tmp_dir, _dir, 'original.prm')
            assert(os.path.isfile(prm_path))
            target_path = os.path.join(results_dir, case_name + '.prm')
            print("copy %s to %s" % (prm_path, target_path))
            shutil.copy(prm_path, target_path)
            statistics_path = os.path.join(results_tmp_dir, _dir, 'statistics')
            assert(os.path.isfile(statistics_path))
            target_path = os.path.join(results_dir, "statistics_" + case_name)
            print("copy %s to %s" % (statistics_path, target_path))
            shutil.copy(statistics_path, target_path)
        else:
            continue

    return results_dir


def analyze_affinity_test_results(test_results_dir, output_dir, **kwargs):
    '''
    analyze affinity test results
    kwargs:
        debug: output debug information
    '''
    debug = kwargs.get("debug", False)
    Statistics = STATISTICS_PLOT('Statistics') # plotter for Statistic outputs
    total_wall_clock = []
    assemble_stokes_system = []
    solve_stokes_system = []
    cores = []
    resolutions = []
    setups = []
    nocs = []
    nofs_stokes = []
    nofs_comp = []
    nofs_temp = []
    # go into sub dirs
    temp_file = os.path.join(ASPECT_LAB_DIR, 'temp')  # file to save partial results
    path_obj = pathlib.Path(test_results_dir).rglob("output*")  # path to all case outputs
    i = 0
    for _path in path_obj:
        i += 1
        output_file = str(_path)
        input_file = os.path.basename(output_file).split("_", maxsplit=1)[1] + ".prm"
        input_path = os.path.join(test_results_dir,  input_file)  # path to the prm file
        statistic_file = "statistics_" + os.path.basename(output_file).split("_", maxsplit=1)[1]
        statistic_path = os.path.join(test_results_dir, statistic_file)
        stokes_solver_type = ParsePrm.GetStokesSolverTypeFromPrm(input_path) 
        patterns = output_file.split('_')
        print("Output file found: %s" % output_file)
        # read and process data
        # 3. Patterns (e.g. name of the server, task per core) are parsed to other arrays (e.g. setups)
        # The data is imported by the parse_block_output.sh
        if stokes_solver_type == "block AMG":
            subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s" 
                            % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        elif stokes_solver_type == "block GMG":
            subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s gmg" 
                            % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        else:
            raise ValueError("stokes_solver_type must be block AMG or block GMG")
        # It is then read by the numpy reader.
        # These numbers are keep in a few arrays (e.g. total_wall_clock).
        try:
            data = np.genfromtxt(temp_file)
            # data dimension:
            # ndim = 1: there is only one step
            # ndim > 1: there are multiple steps
            if data.ndim == 1:
                total_wall_clock.append(data[0])
                assemble_stokes_system.append(data[1])
                solve_stokes_system.append(data[2])
            elif data.ndim == 2:
                total_wall_clock.append(data[0, -1])
                assemble_stokes_system.append(data[1, -1])
                solve_stokes_system.append(data[2, -1])
        except Exception:
            # In the test, this exception means the run fails.
            # It's most like this specific affinity is not
            # good enough for solving this problem at this resolution.
            pass
        else:
            # continue reading other outputs is the wallclock is read successfully
            setups.append(int(patterns[-1]))
            resolutions.append(int(patterns[-2]))
            cores.append(int(patterns[-3]))
            # read data from the statistic files
            Statistics.ReadHeader(statistic_path)
            Statistics.ReadData(statistic_path)
            col_noc = Statistics.header['Number_of_mesh_cells']['col']
            col_dof_stokes = Statistics.header['Number_of_Stokes_degrees_of_freedom']['col']
            col_dof_temperature = Statistics.header['Number_of_temperature_degrees_of_freedom']['col']
            col_dof_composition = Statistics.header['Number_of_degrees_of_freedom_for_all_compositions']['col']
            nocs_raw = Statistics.data[:, col_noc]
            nofs_stokes_raw = Statistics.data[:, col_dof_stokes]
            nofs_temp_raw = Statistics.data[:, col_dof_temperature]
            nofs_comp_raw = Statistics.data[:, col_dof_composition]
            nocs.append(nocs_raw[-1])
            nofs_stokes.append(nofs_stokes_raw[-1])
            nofs_temp.append(nofs_temp_raw[-1])
            nofs_comp.append(nofs_comp_raw[-1])
    # assert there is at least one case output
    Utilities.my_assert(i > 0, AssertionError, "There is no output* file in the folder %s" % test_results_dir)
    # make everything a numpy array
    nocs = np.array(nocs)
    nofs_stokes = np.array(nofs_stokes)
    nofs_temp = np.array(nofs_temp)
    nofs_comp = np.array(nofs_comp)
    setups = np.array(setups)
    resolutions = np.array(resolutions)
    cores = np.array(cores)
    total_wall_clock = np.array(total_wall_clock)
    assemble_stokes_system = np.array(assemble_stokes_system)
    solve_stokes_system = np.array(solve_stokes_system)
    # screen output
    # 1. output numbers from the log.txt file
    # 2. output mesh information from the Statistics file (during loop of resolution)
    if debug:
        print("Affinity Test Results (debug outputs):")
        print("nocs: ", nocs)
        print("nofs_stokes: ", nofs_stokes)
        print("nofs_temp: ", nofs_temp)
        print("nofs_comp: ", nofs_comp)
        print('setups:', setups)
        print('resolutions', resolutions)
        print("cores:", cores)
        print("total wall clocks:", total_wall_clock)
        print("assemble_stokes_system:", assemble_stokes_system)
        print("solve_stokes_system:", solve_stokes_system)
    # plot via matplotlib
    # 1. loop for resolutions, generate one line for every resolution
    # 2. plots are generated in log-log or linear format
    resolution_options = []
    noc_list = []
    nof_stokes_list = []
    nof_temp_list = []
    nof_comp_list = []
    cores_sorted_list = []
    wallclocks_sorted_list = []
    totalclocks_sorted_list = []
    assemble_stokes_system_res_list = []
    solve_stokes_system_res_list = []
    for resolution in resolutions:
        if resolution not in resolution_options:
            resolution_options.append(resolution)
    resolution_options = np.array(resolution_options)
    resolution_options = np.sort(resolution_options)
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        plot_indexes = (resolutions == resolution)
        # parse outputs 
        noc = nocs[plot_indexes][0]
        noc_list.append(float(noc))
        nof_stokes = nofs_stokes[plot_indexes][0]
        nof_stokes_list.append(float(nof_stokes))
        nof_temp = nofs_temp[plot_indexes][0]
        nof_temp_list.append(float(nof_temp))
        nof_comp = nofs_comp[plot_indexes][0]
        nof_comp_list.append(float(nof_comp))
        cores_res = cores[plot_indexes]
        sequential_indexes = np.argsort(cores_res)
        cores_sorted = cores_res[sequential_indexes]
        cores_sorted_list.append(cores_sorted.copy())
        wallclocks_res = total_wall_clock[plot_indexes]
        wallclocks_sorted = wallclocks_res[sequential_indexes]
        wallclocks_sorted_list.append([float(f) for f in wallclocks_sorted])
        totalclocks_sorted_list.append([float(wallclocks_sorted[i]) * float(cores_sorted[i])\
                                        for i in range(len(wallclocks_sorted))])
        assemble_stokes_system_res = assemble_stokes_system[plot_indexes]
        assemble_stokes_system_res_list.append([float(f) for f in assemble_stokes_system_res])
        assemble_stokes_system_sorted = assemble_stokes_system_res[sequential_indexes]
        solve_stokes_system_res = solve_stokes_system[plot_indexes]
        solve_stokes_system_res_list.append([float(f) for f in solve_stokes_system_res])
        solve_stokes_system_sorted = solve_stokes_system_res[sequential_indexes]
        # screen outputs
        print("resolution: %d, cells: %d, degree of freedoms (stokes: %d + temperature: %d + composition: %d) "\
        % (resolution, noc, nof_stokes, nof_temp, nof_comp))
        print("\tnumber of cores: ", cores_sorted)
        # print("\tNumber of Stokes degrees of freedom: ", nofs_stokes[i])
        # print("\tNumber of temperature degrees of freedom: ", nofs_temp[i])
        print("\twallclock (s): ", wallclocks_sorted)
        print("\tassemble stokes system (s): ", assemble_stokes_system_sorted)
        print("\tsolve stokes system (s): ", solve_stokes_system_sorted)
    # plot degree of freedom
    fig, ax= plt.subplots(tight_layout=True, figsize=(5, 5))  # plot of wallclock
    nofs_stokes_list = []
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        label_stokes_dof = "Stokes DOFs(resolution=%d)" % resolution
        label_temperature_dof = "Temperature DOFs(resolution=%d)" % resolution
        label_composition_dof = "Composition DOFs(resolution=%d)" % resolution
        plot_indexes = (resolutions == resolution)
        ax.semilogy(resolution, nofs_stokes[plot_indexes][0], ".",\
                color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_stokes_dof)
        ax.semilogy(resolution, nofs_temp[plot_indexes][0], "o",\
                color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_temperature_dof)
        ax.semilogy(resolution, nofs_comp[plot_indexes][0], "*",\
                color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_composition_dof)
    ax.grid()
    ax.set_xlabel('Cores')
    ax.set_ylabel('DOFs')
    ax.legend()
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_dof.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)
    
    # plot wall clock
    fig = plt.figure(tight_layout=True, figsize=(10, 5))  # plot of wallclock
    gs = gridspec.GridSpec(1, 2)
    ax_wallclock = fig.add_subplot(gs[0, 0])
    ax_wallclock_percentage = fig.add_subplot(gs[0, 1])
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        label_wallclock = 'Total Wall Clock(resolution=%d)' % resolution
        label_wallclock_assemle = 'Assemble Stokes (resolution=%d)' % resolution
        label_wallclock_solve = 'Solve Stokes (resolution=%d)' % resolution
        plot_indexes = (resolutions == resolution)
        cores_res = cores[plot_indexes]
        sequential_indexes = np.argsort(cores_res)
        cores_sorted = cores_res[sequential_indexes]
        wallclocks_res = total_wall_clock[plot_indexes]
        wallclocks_sorted = wallclocks_res[sequential_indexes]
        assemble_stokes_system_res = assemble_stokes_system[plot_indexes]
        assemble_stokes_system_sorted = assemble_stokes_system_res[sequential_indexes]
        solve_stokes_system_res = solve_stokes_system[plot_indexes]
        solve_stokes_system_sorted = solve_stokes_system_res[sequential_indexes]
        ax_wallclock.loglog(cores_sorted, wallclocks_sorted, ".-",\
                            color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_wallclock)
        ax_wallclock.loglog(cores_sorted, assemble_stokes_system_sorted, "--",\
                            color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_wallclock_assemle)
        ax_wallclock.loglog(cores_sorted, solve_stokes_system_sorted, "-.",\
                            color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_wallclock_solve)
        ax_wallclock.set_aspect('equal', adjustable='box')
        ax_wallclock.set_xlabel('Cores')
        ax_wallclock.set_ylabel('Time [s]')
        ax_wallclock.grid()
        ax_wallclock.set_title('Wall Clock')
        ax_wallclock.legend(fontsize='x-small')
        ax_wallclock_percentage.semilogx(cores_sorted, assemble_stokes_system_sorted/wallclocks_sorted, "--",\
                                         color=cm.gist_rainbow(1.0*i/len(resolution_options)))
        ax_wallclock_percentage.semilogx(cores_sorted, solve_stokes_system_sorted/wallclocks_sorted, "-.",\
                                         color=cm.gist_rainbow(1.0*i/len(resolution_options)))
        ax_wallclock_percentage.set_xlabel('Cores')
        ax_wallclock_percentage.set_ylabel('Percentage')
        ax_wallclock_percentage.grid()
        ax_wallclock_percentage.set_title('Percentage of Each Part')
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_wallclock.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)
    # plot wall clock with an ideal line
    fig, ax= plt.subplots(tight_layout=True, figsize=(5, 5))  # plot of wallclock
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        label_wallclock = 'Total Wall Clock(resolution=%d)' % resolution
        plot_indexes = (resolutions == resolution)
        cores_res = cores[plot_indexes]
        sequential_indexes = np.argsort(cores_res)
        cores_sorted = cores_res[sequential_indexes]
        wallclocks_res = total_wall_clock[plot_indexes]
        wallclocks_sorted = wallclocks_res[sequential_indexes]
        ax.loglog(cores_sorted, wallclocks_sorted, ".-",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_wallclock)
        wallclocks_ideal = cores_sorted[0]*wallclocks_sorted[0] / cores_sorted
        ax.loglog(cores_sorted, wallclocks_ideal, "--",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)))
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Cores')
        ax.set_ylabel('Time [s]')
        ax.grid()
        ax.set_title('Wall Clock')
        ax.legend(fontsize='x-small')
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_wallclock_with_ideal.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)
    # plot the speedup
    fig, axs = plt.subplots(1, len(resolution_options), tight_layout=True,\
                            figsize=(5*len(resolution_options), 5))
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        label_speedup = 'Speed Up(resolution=%d)' % resolution
        plot_indexes = (resolutions == resolution)
        cores_res = cores[plot_indexes]
        sequential_indexes = np.argsort(cores_res)
        cores_sorted = cores_res[sequential_indexes]
        wallclocks_res = total_wall_clock[plot_indexes]
        wallclocks_sorted = wallclocks_res[sequential_indexes]
        speedups = wallclocks_sorted[0] / wallclocks_sorted
        speedups_ideal = cores_sorted / cores_sorted[0]
        axs[i].plot(cores_sorted, speedups, ".-",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_speedup)
        axs[i].plot(cores_sorted, speedups_ideal, "--",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)))
    for i in range(len(resolution_options)):
        # add options and lables
        axs[i].set_xlabel('Cores')
        axs[i].set_ylabel('Speed Up')
        axs[i].grid()
        axs[i].legend(fontsize='x-small')
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_speedup.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)

    # plot the parallel efficiency
    fig, axs = plt.subplots(1, len(resolution_options), tight_layout=True,\
                            figsize=(5*len(resolution_options), 5))
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        label_parallel_efficiency = 'Parallel Efficiency(resolution=%d)' % resolution
        plot_indexes = (resolutions == resolution)
        cores_res = cores[plot_indexes]
        sequential_indexes = np.argsort(cores_res)
        cores_sorted = cores_res[sequential_indexes]
        wallclocks_res = total_wall_clock[plot_indexes]
        wallclocks_sorted = wallclocks_res[sequential_indexes]
        parallel_efficiencies = wallclocks_sorted[0]*cores_sorted[0] / wallclocks_sorted / cores_sorted
        parallel_efficiencies_ideal = np.ones(cores_sorted.shape)
        axs[i].plot(cores_sorted, parallel_efficiencies, ".-",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)), label=label_parallel_efficiency)
        axs[i].plot(cores_sorted, parallel_efficiencies_ideal, "--",\
                  color=cm.gist_rainbow(1.0*i/len(resolution_options)))
    for i in range(len(resolution_options)):
        # add options and lables
        axs[i].set_xlabel('Cores')
        axs[i].set_ylabel('Parallel Efficiency')
        axs[i].grid()
        axs[i].legend(fontsize='x-small')
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_parallel_efficiency.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)

    # plot the scale-increasing efficiency
    fig, ax= plt.subplots(tight_layout=True, figsize=(5, 5))  # plot of wallclock
    totalclocks_plist = []
    cores_plist = []
    _label = "cores = ("
    costs = []
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        totalclocks_sorted = np.array(totalclocks_sorted_list[i])
        cores_sorted = cores_sorted_list[i]
        i_min = np.argmin(totalclocks_sorted)
        totalclock_min = totalclocks_sorted[i_min]
        core = cores_sorted[i_min]
        totalclocks_plist.append(totalclock_min)
        cores_plist.append(core)
        if i == 0:
            costs.append(1.0)
        else:
            costs.append(totalclock_min / totalclocks_plist[i-1])
        _label += (" " + str(core))
    _label += ")"
    totalclocks = np.array(totalclocks_plist)
    cores = np.array(cores_plist)
    ax.plot(resolution_options, costs, "o", label=_label)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Refining cost')
    ax.grid()
    ax.legend(fontsize='x-small')
    basename = os.path.basename(test_results_dir)
    fig.tight_layout()
    filepath='%s/%s_refining_cost.png' % (output_dir, basename)
    print("output file generated: ", filepath)
    fig.savefig(filepath)

    # prepare the table of outputs
    # summary table
    tex_contents = ""
    table_contents = ""
    nof_list = [(nof_stokes_list[i] + nof_temp_list[i] + nof_comp_list[i]) for i in range(len(resolution_options))]
    print('nof_list: ', nof_list)  # debug
    header = ["resolution", "number of cells", "dof", "dof stokes", "dof temp", "dof comp"]
    data = [resolution_options, noc_list, nof_list, nof_stokes_list, nof_temp_list, nof_comp_list]
    TexTable = Utilities.TEX_TABLE("table-summary",\
                                    header=header, data=data) # class initiation
    table_contents += TexTable(format="latex", caption="Summary of affinity test")
    tex_contents += (table_contents + '\n')
    # subsequent tables of resolutions
    table_contents = ""
    header = ["number of cores", "wallclock", "assemble stokes system", "solve stokes system", "total clock"]
    for i in range(len(resolution_options)):
        data = []
        resolution = resolution_options[i]
        data.append(cores_sorted_list[i])
        data.append(wallclocks_sorted_list[i])
        data.append(assemble_stokes_system_res_list[i])
        data.append(solve_stokes_system_res_list[i])
        data.append(totalclocks_sorted_list[i])
        TexTable = Utilities.TEX_TABLE("table-%d" % (i+1),\
                                        header=header, data=data) # class initiation
        table_contents += TexTable(format="latex", caption="resolution: %d" % resolution)
        table_contents += '\n'
    tex_contents += table_contents
    filepath='%s/summary.tex' % (output_dir)
    # figures
    _caption = "Wallclock for test %s" % basename # wallclock
    tex_contents += Utilities.latex_figure('%s_wallclock_with_ideal.png' %  basename, caption=_caption)
    tex_contents += "\n\n"
    _caption = "SpeedUp for test %s" % basename # speedup
    tex_contents += Utilities.latex_figure('%s_speedup.png' %  basename, caption=_caption)
    tex_contents += "\n\n"
    _caption = "Parallel Efficiency for test %s" % basename # parallel efficiency
    tex_contents += Utilities.latex_figure('%s_parallel_efficiency.png' %  basename, caption=_caption)
    tex_contents += "\n\n"
    _caption = "Refining cost for test %s" % basename # parallel efficiency
    tex_contents += Utilities.latex_figure('%s_refining_cost.png' %  basename, caption=_caption)
    tex_contents += "\n\n"
    # text
    tex_text_contents = "On date foo, " 
    tex_text_contents += "this test is conducted on the %s partition.\n" % basename
    tex_text_contents += "The model used in this test is foo.\n"
    tex_text_contents += "These cases are run for foo steps.\n"
    tex_text_contents += "An parameter file and a slurm file of this case could be found in folder foo.\n"
    tex_text_contents += "The original path of this test is %s.\n" % test_results_dir
    tex_text_contents += "wallclock (figure \\ref{fig:%s_wallclock_with_ideal});\n" % basename
    tex_text_contents += "speedup (figure \\ref{fig:%s_speedup});\n" % basename
    tex_text_contents += "parallel efficiency (figure \\ref{fig:%s_parallel_efficiency});\n" % basename
    tex_text_contents += "refining cost (figure \\ref{fig:%s_refining_cost});\n" % basename
    tex_text_contents += "Based on these results, we suggest running the model on foo nodes foo cores.\n"
    tex_text_contents += "The resources needed for one case is estimated to be foo.\n"
    tex_contents = tex_text_contents + "\n" + tex_contents
    with open(filepath, 'w')  as fileout:
        fileout.write(tex_contents)
    print("latex summary generated: ", filepath)


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
                        project=Affinity_Opt.get_project(), build_directory=Affinity_Opt.get_build_directory(),\
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
        create_tests_with_json(arg.json, AFFINITY, AFFINITY_OPT)

    elif commend == 'analyze_results':
        # example:
        # python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        # -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/affinity_test_20211025 -c peloton-rome-128tasks-socket-openmpi-4.1.0
        assert(os.path.isdir(arg.inputs))  # check output dir exists
        test_results_dir = organize_result(arg.inputs, arg.cluster)
        img_dir = os.path.join(arg.inputs, "img")
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        analyze_affinity_test_results(test_results_dir, img_dir)
    elif commend in ['help', '-h']:
        Usage()
    else:
        raise ValueError("No such command (%s), see Usage (run with help or -h)" % commend)

# run script
if __name__ == '__main__':
    main()
