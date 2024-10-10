# -*- coding: utf-8 -*-
r"""Plot visit output

This exports:

  -

This depends on:

  -
Examples of usage:

  - default usage:

        python -m

descriptions

Notes:
    on the name standard:
        step - computational step in model
        vtu_step - indexing of steps of vtu outputs. (e.g. step = 0, vtu_step = 0; step = 1, vtu_step = 100)
        vtu_snapshot - indexing of output files of vtu (e.g. solution-00001.pvtu, then vtu_snapshot = 1)
"""
import numpy as np
import sys, os, argparse
import json, re
# import pathlib
import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Plot as Plot
import shilofue.ParsePrm as ParsePrm
from shilofue.CaseOptions import CASE_OPTIONS
from shilofue.PlotDepthAverage import ExportData
import warnings

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
Visualization with visit\n\
\n\
Examples of usage: \n\
\n\
  - translate script: \n\
\n\
        Lib_PlotVisit visit_options -i $TwoDSubduction_DIR/latent_heat_issue/cookbook_latent-heat -sr temperature.py\n\
            -sr: the relative path under the visit_script folder\n\
\n\
  - run script: \n\
        Lib_PlotVisit run -i $TwoDSubduction_DIR/non_linear34/eba_low_tol_newton_shift_CFL0.8/visit_scripts/slab.py\n\
\n\
  - run vtk scripts: \n\
        Lib_PlotVisit vtk_options -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8_lh -p TwoDSubduction_MOW -s 25\n\
            -p operations, available operations are: \n\
                TwoDSubduction_MOW: pull out tentative MOW from temperature: \n\
                TwoDSubduction_SlabAnalysis: analyze slab morphology: \n\
        ")


class VISIT_OPTIONS(CASE_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
            steps (int): plot some steps
            last_step(list): plot the last few steps
        """
        steps = kwargs.get('steps', None)
        last_step = kwargs.get('last_step', None)
        time_interval = kwargs.get('time_interval', None)
        plot_axis = kwargs.get('plot_axis', False)
        max_velocity = kwargs.get('max_velocity', -1.0)
        slices = kwargs.get('slices', 3)
        graphical_type = kwargs.get("graphical_type", "pvd")
        # call function from parent
        CASE_OPTIONS.Interpret(self)
        # particle file
        particle_file = os.path.join(self._output_dir, 'particles.visit')
        if os.access(particle_file, os.R_OK):
            self.options["VISIT_PARTICLE_FILE"] = particle_file
        # visit file
        self.options["VISIT_FILE"] = self.visit_file
        self.options["PARAVIEW_FILE"] = self.paraview_file
        # data types
        self.options["HAS_DYNAMIC_PRESSURE"] = '0'
        try:
            visualization_output_variables = self.idict['Postprocess']['Visualization']['List of output variables']
        except KeyError:
            pass
        else:
            if re.match('.*nonadiabatic\ pressure', visualization_output_variables):
                self.options["HAS_DYNAMIC_PRESSURE"] = '1'

        # plot options
        # plot axis
        if plot_axis:
            self.options["PLOT_AXIS"] = '1'
        else: 
            self.options["PLOT_AXIS"] = '0'
        # maximum velocity
        self.options["MAX_VELOCITY"] = str(max_velocity)
        self.options["PLOT_TYPES"] = str(kwargs.get('plot_types', []))
        # additional fields to load for model
        additional_fields = kwargs.get('additional_fields', [])
        assert(type(additional_fields) == list)
        self.options["ADDITIONAL_FIELDS"] = str(additional_fields)

        # get all the available snaps for ploting by checking on the existence of the pvtu file
        # the correspondent, time, time step are also figured out.
        graphical_snaps_guess, times_guess, time_steps_guess = GetSnapsSteps(self._case_dir, 'graphical')
        graphical_snaps = []
        time_steps = []
        times = []
        for i in range(len(graphical_snaps_guess)):
            graphical_snap = graphical_snaps_guess[i]
            time_step = time_steps_guess[i]
            _time = times_guess[i]
            graphical_file_path = None
            if graphical_type == "pvd":
                graphical_file_path = os.path.join(self.options["DATA_OUTPUT_DIR"], "solution", "solution-%05d.pvtu" % graphical_snap)
            elif graphical_type == "slice_center":
                graphical_file_path = os.path.join(self._case_dir, "vtk_outputs", "center_profile_%05d.txt" % graphical_snap)
            if os.path.isfile(graphical_file_path):
                graphical_snaps.append(graphical_snap)
                time_steps.append(time_step)
                times.append(_time)
        self.all_graphical_snaps = graphical_snaps
        self.all_graphical_timesteps = time_steps 
        self.all_graphical_times = times
        self.options['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'] = str(graphical_snaps)
        self.options['ALL_AVAILABLE_GRAPHICAL_TIMESTEPS'] = str(time_steps)
        self.options['ALL_AVAILABLE_GRAPHICAL_TIMES'] = str(times)
        particle_snaps, _, _ = GetSnapsSteps(self._case_dir, 'particle')
        self.options['ALL_AVAILABLE_PARTICLE_SNAPSHOTS'] = str(particle_snaps)
        particle_output_dir = os.path.join(self._output_dir, "slab_morphs")
        self.options["PARTICLE_OUTPUT_DIR"] = particle_output_dir
        # get the last step in the series
        try:
            self.last_step = max(0, graphical_snaps[-1] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT']))  # it is the last step we have outputs
        except IndexError:
            # no snaps, stay on the safe side
            self.last_step = -1
        # add an option of the last step
        self.options["LAST_STEP"] = self.last_step

        # set steps to plot
        # Priority:
        #   1. a list of steps
        #   2. the last few steps
        #   3. only the last step
        if type(steps) == list:
            for step in steps:
                assert(type(step) == int)
            self.options['GRAPHICAL_STEPS'] = steps  # always plot the 0 th step
        elif type(time_interval) is float and time_interval > 0.0:
            times_ndarray = np.array(times)
            time_series_from_interval = np.arange(0.0, times[-1], time_interval, dtype=float)
            self.options['GRAPHICAL_STEPS'] = []
            for i in range(time_series_from_interval.size):
                _time = time_series_from_interval[i]
                idx = np.argmin(abs(times_ndarray - _time))
                self.options['GRAPHICAL_STEPS'].append(graphical_snaps[idx] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT']))
        elif type(last_step) == int:
            # by this option, plot the last few steps
            self.options['GRAPHICAL_STEPS'] = [0]  # always plot the 0 th step
            self.options['GRAPHICAL_STEPS'] += [i for i in range(max(self.last_step - last_step + 1, 0), self.last_step + 1)]
        elif type(steps) == str and steps == "auto":
            # 
            # determine the options by the number of steps and slice them by the number of slices
            assert(slices > 0)
            self.options['GRAPHICAL_STEPS'] = [int(i) for i in np.linspace(0 , int(self.options["LAST_STEP"]), slices)]
            # self.options['GRAPHICAL_STEPS'].append(int(self.options["LAST_STEP"]))
        else:
            # by default append the first and the computing step.
            self.options['GRAPHICAL_STEPS'] = [0]
            if self.last_step > 0:
                self.options['GRAPHICAL_STEPS'].append(self.last_step)

        # get time steps
        self.options['GRAPHICAL_TIME_STEPS'] = []
        for step in self.options['GRAPHICAL_STEPS']:
            found = False
            for i in range(len(graphical_snaps)):
                if step == max(0, graphical_snaps[i] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])):
                    found = True
                    self.options['GRAPHICAL_TIME_STEPS'].append(time_steps[i])
            if not found:
                warnings.warn("%s: step %d is not found" % (Utilities.func_name(), step))
            # Utilities.my_assert(found, ValueError, "%s: step %d is not found" % (Utilities.func_name(), step))

    def visit_options(self, extra_options):
        '''
        deprecated
        '''
        # optional settings
        for key, value in extra_options.items():
            # slab
            if key == 'slab':
                self.options['IF_PLOT_SLAB'] = 'True'
                self.options['GRAPHICAL_STEPS'] = value.get('steps', [0])
                self.options['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
            # export particles for slab morph
            elif key == 'slab_morph':
                self.options['IF_EXPORT_SLAB_MORPH'] = 'True'
                # check directory
                if not os.path.isdir(particle_output_dir):
                    os.mkdir(particle_output_dir)
    
    def vtk_options(self, **kwargs):
        '''
        options of vtk scripts
        '''
        generate_horiz_file = kwargs.get('generate_horiz', False)
        operation = kwargs.get('operation', 'default')
        vtu_step = int(kwargs.get('vtu_step', 0))
        # houriz_avg file
        if generate_horiz_file:
            _time, time_step = self.get_time_and_step(vtu_step)
            depth_average_path = os.path.join(self.options["DATA_OUTPUT_DIR"], 'depth_average.txt')
            assert(os.path.isfile(depth_average_path))
            output_dir = os.path.join(self._case_dir, 'temp_output')
            try:  # This works better in parallel
                os.mkdir(output_dir)
            except FileExistsError:
                pass
            _, ha_output_file = ExportData(depth_average_path, output_dir, time_step=time_step, fix_time_step=True)
            self.options['VTK_HORIZ_FILE'] = ha_output_file
        else:
            self.options['VTK_HORIZ_FILE'] = os.path.join(ASPECT_LAB_DIR, 'output', 'depth_average_output')
        # directory to output from vtk script
        self.options['VTK_OUTPUT_DIR'] = os.path.join(self._case_dir, "vtk_outputs")
        if not os.path.isdir(self.options['VTK_OUTPUT_DIR']):
            os.mkdir(self.options['VTK_OUTPUT_DIR'])
        # file to read in vtk
        Utilities.my_assert((vtu_step >= 0 and vtu_step <= self.last_step), ValueError, "vtu_step needs to be within the range of [%d, %d]" % (0, self.last_step))  # check the range of steps
        self.options['PVTU_FILE'] = os.path.join(self._output_dir, "solution", "solution-%05d.pvtu" % (vtu_step + int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])))
        # type of operation
        self.options['OPERATION'] = operation

    def get_time_and_step(self, vtu_step):
        '''
        Convert vtu_step to step and time in model
        ''' 
        assert(len(self.all_graphical_snaps) > 0)
        assert(len(self.all_graphical_timesteps) > 0)
        # find step in all available steps
        found = False
        i = 0
        for snap_shot in self.all_graphical_snaps:
            if vtu_step == max(0, int(snap_shot) - int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])):
                found = True
                step = int(self.all_graphical_timesteps[i])
            i += 1
        Utilities.my_assert(found, ValueError, "%s: vtu_step %d is not found" % (Utilities.func_name(), vtu_step))
        time = self.Statistics.GetTime(step)
        return time, step
    
    def get_time_and_step_by_snapshot(self, vtu_snapshot):
        '''
        Convert vtu_snapshot to step and time in model
        ''' 
        assert(len(self.all_graphical_snaps) > 0)
        assert(len(self.all_graphical_timesteps) > 0)
        # find step in all available steps
        found = False
        i = 0
        for snap_shot in self.all_graphical_snaps:
            if vtu_snapshot == snap_shot:
                found = True
                step = int(self.all_graphical_timesteps[i])
            i += 1
        Utilities.my_assert(found, ValueError, "%s: vtu_snapshot %d is not found" % (Utilities.func_name(), vtu_snapshot))
        time = self.Statistics.GetTime(step)
        return time, step

    def get_timestep_by_time(self, _time: float):
        '''
        Retrieves the closest graphical time and its corresponding timestep based on the given time.
        
        Parameters:
        _time (float): The reference time for which the closest graphical time and timestep are to be found.
        
        Returns:
        Tuple: A tuple containing the closest graphical time and its associated timestep.
        '''
        index = np.argmin(np.abs(np.array(self.all_graphical_times) - _time))
        return self.all_graphical_times[index], self.all_graphical_timesteps[index], self.all_graphical_snaps[index] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])



class PARALLEL_WRAPPER_FOR_VTK():
    '''
    a parallel wrapper for analyzing slab morphology
    Attributes:
        name(str): name of this plot
        case_dir (str): case directory
        module (a function): a function to use for plotting
        last_pvtu_step (str): restart from this step, as there was previous results
        if_rewrite (True or False): rewrite previous results if this is true
        pvtu_steps (list of int): record the steps
        outputs (list of str): outputs
    '''
    def __init__(self, name, module, **kwargs):
        '''
        Initiation
        Inputs:
            name(str): name of this plot
            module (a function): a function to use for plotting
            kwargs (dict)
                last_pvtu_step
                if_rewrite
        '''
        self.name = name
        self.module = module
        self.last_pvtu_step = kwargs.get('last_pvtu_step', -1)
        self.if_rewrite = kwargs.get('if_rewrite', False)
        self.do_assemble = kwargs.get('assemble', True)
        self.kwargs = kwargs
        self.pvtu_steps = []
        self.outputs = []
        pass

    def configure(self, case_dir):
        '''
        configure
        Inputs:
            case_dir (str): case diretory to assign
        '''
        os.path.isdir(case_dir)
        self.case_dir = case_dir
    
    def __call__(self, pvtu_step):
        '''
        call function
        Inputs:
            pvtu_step (int): the step to plot
        '''
        expect_result_file = os.path.join(self.case_dir, 'vtk_outputs', '%s_s%06d' % (self.name, pvtu_step))
        if pvtu_step <= self.last_pvtu_step and not self.if_rewrite:
            # skip existing steps
            return 0
        if os.path.isfile(expect_result_file) and not self.if_rewrite:
            # load file content
            print("%s: previous result exists(%s), load" % (Utilities.func_name(), expect_result_file))
            with open(expect_result_file, 'r') as fin:
                pvtu_step = int(fin.readline())
                output = fin.readline()
        else:
            if self.do_assemble:    
                # here the outputs from individual steps are combined together
                if pvtu_step == 0:
                    # start new file with the 0th step
                    pvtu_step, output = self.module(self.case_dir, pvtu_step, new=True, **self.kwargs)
                else:
                    pvtu_step, output = self.module(self.case_dir, pvtu_step, **self.kwargs)
                with open(expect_result_file, 'w') as fout:
                    fout.write('%d\n' % pvtu_step)
                    fout.write(output)
                print("%s: pvtu_step - %d, output - %s" % (Utilities.func_name(), pvtu_step, output))
                # self.pvtu_steps.append(pvtu_step) # append to data
                self.outputs.append(output)
            else: 
                # otherwise, just call the module for each steps
                if pvtu_step == 0:
                    # start new file with the 0th step
                    self.module(self.case_dir, pvtu_step, new=True, **self.kwargs)
                else:
                    self.module(self.case_dir, pvtu_step, **self.kwargs)
        return 0
    
    def assemble(self):
        '''
        Returns:
            pvtu_steps
            outputs
        '''
        assert(len(self.pvtu_steps) == len(self.outputs))
        length = len(self.pvtu_steps)
        # bubble sort
        for i in range(length):
            for j in range(i+1, length):
                if self.pvtu_steps[j] < self.pvtu_steps[i]:
                    temp = self.pvtu_steps[i]
                    self.pvtu_steps[i] = self.pvtu_steps[j]
                    self.pvtu_steps[j] = temp
                    temp = self.outputs[i]
                    self.outputs[i] = self.outputs[j]
                    self.outputs[j] = temp
        return self.pvtu_steps, self.outputs
    
    def assemble_parallel(self):
        '''
        Returns:
            pvtu_steps
            outputs
        '''
        for pvtu_step in self.pvtu_steps:
            expect_result_file = os.path.join(self.case_dir, 'vtk_outputs', '%s_s%06d' % (self.name, pvtu_step))
            assert(os.path.isfile(expect_result_file))
            with open(expect_result_file, 'r') as fin:
                fin.readline()
                output = fin.readline()
                self.outputs.append(output)
        return self.pvtu_steps, self.outputs
    
    def set_pvtu_steps(self, pvtu_steps):
        '''
        set_pvtu_steps
        Inputs:
            pvtu_steps(list of int): step to look for
        '''
        self.pvtu_steps = pvtu_steps
    
    def delete_temp_files(self, pvtu_steps):
        '''
        delete temp files
        Inputs:
            pvtu_steps(list of int): step to look for
        '''
        print('delete temp files')
        for pvtu_step in pvtu_steps:
            expect_result_file = os.path.join(self.case_dir, 'vtk_outputs', '%s_s%06d' % (self.name, pvtu_step))
            if os.path.isfile(expect_result_file):
                os.remove(expect_result_file)
    
    def clear(self):
        '''
        clear data
        '''
        self.pvtu_steps = []
        self.outputs = []


class PREPARE_RESULT_OPTIONS(CASE_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
        """
        # call function from parent
        CASE_OPTIONS.Interpret(self)
        visual_step = kwargs.get('step', 0)
        # step & format for a visit output
        self.options['NUMERICAL_STEP'] = "%06d" % visual_step
        self.options['GRAPHICAL_STEP'] =  "%06d" % (visual_step + int(self.options['INITIAL_ADAPTIVE_REFINEMENT']))
        # time
        try:
            time_between_graphical_output = float(self.idict['Postprocess']['Visualization']['Time between graphical output'])
        except KeyError:
            time_between_graphical_output = 1e8
        time = visual_step * time_between_graphical_output
        step = self.Statistics.GetStep(time)
        self.options['STEP_TIME_STAMP'] = "step %d, %.4e yrs" % (step, time)
        self.options['NUMERICAL_TIME_4E'] = "%.4e" % time

 
def GetSnapsSteps(case_dir, type_='graphical'):
    '''
    Get snaps for visualization from the record of statistic file.
    This function requires a consistent statistic file with respect to the vtu files.
    Checking the statistic file is more on the safe side from checking the vtu file, since a newer run might overight the older result.
    '''
    case_output_dir = os.path.join(case_dir, 'output')

    # import parameters
    prm_file = os.path.join(case_dir, 'output', 'original.prm')
    Utilities.my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,
              'case prm file - %s cannot be read' % prm_file)
    with open(prm_file, 'r') as fin:
        idict = ParsePrm.ParseFromDealiiInput(fin)

    # import statistics file
    Statistics = Plot.STATISTICS_PLOT_OLD('Statistics')
    statistic_file = os.path.join(case_output_dir, 'statistics')
    Utilities.my_assert(os.access(statistic_file, os.R_OK), FileNotFoundError,
              'case statistic file - %s cannot be read' % prm_file)
    Statistics.ReadHeader(statistic_file)
    Statistics.ReadData(statistic_file)
    col_time = Statistics.header['Time']['col']
    col_step = Statistics.header['Time_step_number']['col']

    # final time and step
    final_time = Statistics.data[-1, col_time]
    final_step = int(Statistics.data[-1, col_step])

    total_graphical_outputs = 0
    graphical_times = []
    graphical_steps = []
    # time interval
    # graphical
    try:
        time_between_graphical_output = float(idict['Postprocess']['Visualization']['Time between graphical output'])
    except KeyError:
        time_between_graphical_output = 1e8
    if time_between_graphical_output < 1e-6:
        # in case of 0, results are written every step
        total_graphical_outputs = int(final_step) + 1
        graphical_times = Statistics.data[:, col_time]
        graphical_steps = range(final_step + 1)
    else:
        total_graphical_outputs = int(final_time / time_between_graphical_output) + 1
        graphical_times = [i*time_between_graphical_output for i in range(total_graphical_outputs)]
        graphical_steps = [Statistics.GetStep(time) for time in graphical_times]
    # particle
    try:
        time_between_particles_output = float(idict['Postprocess']['Particles']['Time between data output'])
        total_particles_outputs = int(final_time / time_between_particles_output) + 1
    except KeyError:
        time_between_particles_output = 1e8
        total_particles_outputs = 0
    particle_times = [i*time_between_particles_output for i in range(total_particles_outputs)]
    particle_steps = [Statistics.GetStep(time) for time in particle_times]

    # initial_snap
    try:
        initial_snap = int(idict['Mesh refinement']['Initial adaptive refinement'])
    except KeyError:
        initial_snap = 0

    # end snap
    snaps = [0]
    if type_ == 'graphical':
        start_ = initial_snap
        end_ = total_graphical_outputs + initial_snap
        snaps = list(range(start_, end_))
        times = graphical_times
        steps = graphical_steps
    elif type_ == 'particle':
        start_ = 0
        end_ = total_particles_outputs
        snaps = list(range(start_, end_))
        times = particle_times
        steps = particle_steps

    return snaps, times, steps


def RunScripts(visit_script):
    '''
    run script in visit
    Inputs:
    '''
    os.system("echo \"exit()\" | eval \"visit -nowin -cli -s %s\"" % visit_script)


def PrepareVTKOptions(VISIT_OPTIONS, case_dir, type, **kwargs):
    '''
    prepare vtk options for vtk scripts
    Inputs:
        VISIT_OPTIONS: class for the options of visit
        kwargs:
            output
            vtu_step
            include_step_in_filename
    '''
    vtu_step = kwargs.get('vtu_step', 0)
    generate_horiz = kwargs.get('generate_horiz', False)
    operation = kwargs.get('operation', 'default')
    include_step_in_filename = kwargs.get('include_step_in_filename', False)
    vtk_config_dir = os.path.join(ASPECT_LAB_DIR, 'vtk_scripts', "inputs")
    assert(os.path.isdir(vtk_config_dir))
    vtk_config_file = os.path.join(vtk_config_dir, "%s.input" % type)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    Visit_Options.vtk_options(vtu_step=vtu_step, generate_horiz=generate_horiz, operation=operation)
    Visit_Options.read_contents(vtk_config_file)
    Visit_Options.substitute()
    try:
        ofile = kwargs['output']
    except KeyError:
        if include_step_in_filename:
            ofile = os.path.join(Visit_Options.options['VTK_OUTPUT_DIR'], "%s_s%06d" % (os.path.basename(vtk_config_file), vtu_step))
        else:
            ofile = os.path.join(Visit_Options.options['VTK_OUTPUT_DIR'], os.path.basename(vtk_config_file))
    ofile_path = Visit_Options.save(ofile)
    print('%s: %s generated' % (Utilities.func_name(), ofile_path))
    # get time and step
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    return ofile_path, _time, step


def RunVTKScripts(operation, vtk_option_path):
    '''
    run script of vtk
    Inputs:
    Return:
        filename (str): file generated from a vtk script
    '''
    vtk_executable = os.path.join(ASPECT_LAB_DIR, 'vtk_scripts', 'build', operation)
    print("run vtk script: \n%s %s" % (vtk_executable, vtk_option_path))
    completed_process = subprocess.run([vtk_executable, vtk_option_path], capture_output=True, text=True)
    print("\nOutput from vtk script:\n %s" % completed_process.stdout)
    return completed_process.stdout


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    parser.add_argument('-sr', '--script', type=str,
                        default='temperature.py',
                        help='script')
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-p', '--operation', type=str,
                        default='',
                        help='operation to take')
    parser.add_argument('-s', '--step', type=int,
                        default='0',
                        help='step')
    parser.add_argument('-ls', '--last_step', type=bool,
                        default=False,
                        help='Only plot last step (bool value)')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend in ['-h', '--help']:
        # example:
        Usage()

    elif _commend == 'visit_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        case_dir = arg.inputs
        Visit_Options = VISIT_OPTIONS(case_dir)
        # call function
        Visit_Options.Interpret()
        # ofile = os.path.join('visit_scripts', 'slab_sph.py')
        ofile = os.path.join('visit_scripts', os.path.basename(arg.script))
        visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', arg.script)
        visit_base_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'base.py')  # base.py : base file
        Visit_Options.read_contents(visit_base_script, visit_script)  # this part combines two scripts
        Visit_Options.substitute()  # substitute keys in these combined file with values determined by Interpret() function
        ofile_path = Visit_Options.save(ofile, relative=True)  # save the altered script
        pass
    
    elif _commend == 'vtk_options':
        vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, arg.inputs, arg.operation, step=arg.step)
        RunVTKScripts(arg.operation, vtk_option_path)
    
    elif _commend == 'run':
        RunScripts(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
