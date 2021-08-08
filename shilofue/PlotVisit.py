# -*- coding: utf-8 -*-
r"""(one line description)

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
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Utilities as Utilities
import shilofue.ParsePrm as ParsePrm
import shilofue.Plot as Plot

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \
        ")


class VISIT_OPTIONS(ParsePrm.CASE_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, kwargs={}):
        """
        Interpret the inputs, to be reloaded in children
        """
        # call function from parent
        ParsePrm.CASE_OPTIONS.Interpret(self, kwargs)

        # visit file
        self.odict["VISIT_FILE"] = self._visit_file

        # particle file
        particle_file = os.path.join(self._output_dir, 'particles.visit')
        if os.access(particle_file, os.R_OK):
            self.odict["VISIT_PARTICLE_FILE"] = particle_file

        # directory to output data
        self.odict["DATA_OUTPUT_DIR"] = self._output_dir

        # directory to output images
        if not os.path.isdir(self._img_dir):
            os.mkdir(self._img_dir)
        self.odict["IMG_OUTPUT_DIR"] = self._img_dir

        # own implementations
        # initial adaptive refinement
        self.odict['INITIAL_ADAPTIVE_REFINEMENT'] = self.idict['Mesh refinement'].get('Initial adaptive refinement', '6')

        # get snaps for plots
        graphical_snaps, _, _ = GetSnapsSteps(self._case_dir, 'graphical')
        self.odict['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'] = str(graphical_snaps)
        particle_snaps, _, _ = GetSnapsSteps(self._case_dir, 'particle')
        self.odict['ALL_AVAILABLE_PARTICLE_SNAPSHOTS'] = str(particle_snaps)


class PLOT_VISIT(VISIT_OPTIONS):
    """
    inherite the VISIT_OPTIONS clase from Parse.py
    in order to add in additional settings
    """
    # todo
    def Interpret(self, kwargs={}):
        """
        Interpret the inputs, to be reloaded in children
        """
        # call function from parent
        VISIT_OPTIONS.Interpret(self, kwargs)

        # default settings
        self.odict['IF_PLOT_SLAB'] = 'False'
        self.odict['IF_EXPORT_SLAB_MORPH'] = 'False'
        particle_output_dir = os.path.join(self._output_dir, "slab_morphs")
        self.odict["PARTICLE_OUTPUT_DIR"] = particle_output_dir

        # optional settings
        for key, value in kwargs.items():
            # slab
            if key == 'slab':
                self.odict['IF_PLOT_SLAB'] = 'True'
                self.odict['PLOT_SLAB_STEPS'] = value.get('steps', [0])
                self.odict['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
            # export particles for slab morph
            elif key == 'slab_morph':
                self.odict['IF_EXPORT_SLAB_MORPH'] = 'True'
                # check directory
                if not os.path.isdir(particle_output_dir):
                    os.mkdir(particle_output_dir)


def GetSnapsSteps(case_dir, type_='graphical'):
    case_output_dir = os.path.join(case_dir, 'output')

    # import parameters
    prm_file = os.path.join(case_dir, 'case.prm')
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
    except KeyError:
        time_between_particles_output = 1e8
    total_particles_outputs = int(final_time / time_between_particles_output) + 1
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


def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
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
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()

    elif _commend == 'visit_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        # todo
        case_dir = arg.inputs

        Visit_Options = PLOT_VISIT(case_dir)

        # load extra options
        if arg.json_file == './config_case.json':
            # no json file is giving
            extra_options = {}
        else:
            with open(arg.json_file, 'r') as fin:
                dict_in = json.load(fin)
                extra_options = dict_in.get('visit', {})

        # call function
        ofile = os.path.join(ASPECT_LAB_DIR, 'visit_keys_values')
        Visit_Options(ofile, extra_options)
        pass

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
