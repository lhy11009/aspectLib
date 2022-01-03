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
        # call function from parent
        CASE_OPTIONS.Interpret(self)
        # particle file
        particle_file = os.path.join(self._output_dir, 'particles.visit')
        if os.access(particle_file, os.R_OK):
            self.options["VISIT_PARTICLE_FILE"] = particle_file
        # visit file
        self.options["VISIT_FILE"] = self._visit_file
        # houriz_avg file
        self.options['VTK_HORIZ_FILE'] = os.path.join(ASPECT_LAB_DIR, 'output', 'depth_average_output')
        # get snaps for plots
        graphical_snaps_guess, _, _ = GetSnapsSteps(self._case_dir, 'graphical')
        graphical_snaps = []
        for graphical_snap in graphical_snaps_guess:
            pvtu_file_path = os.path.join(self.options["DATA_OUTPUT_DIR"], "solution", "solution-%05d.pvtu" % graphical_snap)
            if os.path.isfile(pvtu_file_path):
                graphical_snaps.append(graphical_snap)
        
        self.options['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'] = str(graphical_snaps)
        particle_snaps, _, _ = GetSnapsSteps(self._case_dir, 'particle')
        self.options['ALL_AVAILABLE_PARTICLE_SNAPSHOTS'] = str(particle_snaps)
        
        particle_output_dir = os.path.join(self._output_dir, "slab_morphs")
        self.options["PARTICLE_OUTPUT_DIR"] = particle_output_dir
        try:
            self.last_step = graphical_snaps[-1] - int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])  # it is the last step we have outputs
        except IndexError:
            # no snaps, stay on the safe side
            self.last_step = -1
        steps = kwargs.get('steps', None)
        last_step = kwargs.get('last_step', None)
        # set steps to plot
        if type(steps) == list:
            for step in steps:
                assert(type(step) == int)
            self.options['GRAPHICAL_STEPS'] = steps  # always plot the 0 th step
        elif type(last_step) == int:
            # by this option, plot the last few steps
            self.options['GRAPHICAL_STEPS'] = [0]  # always plot the 0 th step
            self.options['GRAPHICAL_STEPS'] += [i for i in range(self.last_step - last_step + 1, self.last_step + 1)]
        else:
            self.options['GRAPHICAL_STEPS'] = [0, 1, 2, 3, 4, 5, 6, 7]

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
        operation = kwargs.get('operation', 'slab')
        vtk_step = int(kwargs.get('vtk_step', 0))
        # directory to output from vtk script
        self.options['VTK_OUTPUT_DIR'] = os.path.join(self._case_dir, "vtk_outputs")
        if not os.path.isdir(self.options['VTK_OUTPUT_DIR']):
            os.mkdir(self.options['VTK_OUTPUT_DIR'])
        # file to read in vtk
        Utilities.my_assert((vtk_step >= 0 and vtk_step <= self.last_step), ValueError, "vtk_step needs to be within the range of [%d, %d]" % (0, self.last_step))  # check the range of steps
        self.options['PVTU_FILE'] = os.path.join(self._output_dir, "solution", "solution-%05d.pvtu" % (vtk_step + int(self.options['INITIAL_ADAPTIVE_REFINEMENT'])))

    def get_time_and_step(self, vtk_step):
        '''
        Convert vtk_step to step and time in model
        ''' 
        try:
            time_between_graphical_output = float(self.idict['Postprocess']['Visualization']['Time between graphical output'])
        except KeyError:
            time_between_graphical_output = 1e8
        _time = vtk_step * time_between_graphical_output
        step = self.Statistics.GetStep(_time)
        return _time, step


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
        step = kwargs.get('step', 0)
        # step & format for a visit output
        self.options['NUMERICAL_STEP'] = "%06d" % step
        self.options['GRAPHICAL_STEP'] =  "%06d" % (step + int(self.options['INITIAL_ADAPTIVE_REFINEMENT']))

 
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

    # final time
    final_time = Statistics.data[-1, col_time]

    # time interval
    # graphical
    try:
        time_between_graphical_output = float(idict['Postprocess']['Visualization']['Time between graphical output'])
    except KeyError:
        time_between_graphical_output = 1e8
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
        initial_snap = 6

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


def PrepareVTKOptions(case_dir, operation, **kwargs):
    '''
    prepare vtk options for vtk scripts
    '''
    vtk_config_dir = os.path.join(ASPECT_LAB_DIR, 'vtk_scripts', "inputs")
    assert(os.path.isdir(vtk_config_dir))
    vtk_config_file = os.path.join(vtk_config_dir, "%s.input" % operation)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    vtk_step = kwargs.get('vtk_step', 0)
    Visit_Options.vtk_options(vtk_step=vtk_step)
    Visit_Options.read_contents(vtk_config_file)
    Visit_Options.substitute()
    try:
        ofile = kwargs['output']
    except KeyError:
        ofile = os.path.join(Visit_Options.options['VTK_OUTPUT_DIR'], os.path.basename(vtk_config_file))
    ofile_path = Visit_Options.save(ofile)
    print('%s: %s generated' % (Utilities.func_name(), ofile_path))
    # get time and step
    _time, step = Visit_Options.get_time_and_step(vtk_step)
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
        vtk_option_path, _, _ = PrepareVTKOptions(arg.inputs, arg.operation, step=arg.step)
        RunVTKScripts(arg.operation, vtk_option_path)
    
    elif _commend == 'run':
        RunScripts(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
