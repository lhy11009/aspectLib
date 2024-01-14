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
# import json, re
# import pathlib
# import subprocess
import shilofue.Group as GroupP
from shilofue.Group import CreateGroup
from shilofue.TwoDSubduction0.Cases import CASE, CASE_OPT
from shilofue.TwoDSubduction0.VtkPp import SLABPLOT
from shilofue.Plot import LINEARPLOT
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage, create group: \n\
        Lib_TwoDSubduction0_Group create_group -j ~/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/TwoDSubduction/test_group/test.json\n"
        )

def ShowJsonOption():
    Group_Opt = GROUP_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Group_Opt.document_str()
        )


class GROUP_OPT(GroupP.GROUP_OPT):
    pass

class GROUP(GroupP.GROUP):
    pass


class CASE_SUMMARY(GroupP.CASE_SUMMARY):
    '''
    Attributes:
        cases: name of cases
        steps: end steps of cases
        times: end times of cases
        wallclocks: running time of cases on the wall clock
        ab_paths: absolution_paths of cases
        t660s: time the slab tip reaches 660 km
    '''
    def __init__(self, **kwargs):
        '''
        initiation
        Inputs:
            kwargs
        '''
        GroupP.CASE_SUMMARY.__init__(self, **kwargs)
        self.t660s = []
        self.attrs.append("t660s")
        self.sz_thicks = []
        self.attrs.append("sz_thicks")
        self.sz_depths = []
        self.attrs.append("sz_depths")
        self.sz_viscs = []
        self.attrs.append("sz_viscs")
        self.if_subducts = []
        self.attrs.append("if_subducts")
    
    def Update(self, **kwargs):
        '''
        Update on properties
        Inputs:
            kwargs:
                actions (list): actions to take
        '''
        actions = kwargs.get('actions', [])

        GroupP.CASE_SUMMARY.Update(self, **kwargs)
        
        if "t660" in actions:
            # initiate these field
            self.t660s = [-1 for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_t660(i)

        if "shear_zone" in actions:
            # initiate these field
            self.sz_thicks = [-1 for i in range(self.n_case)]
            self.sz_depths = [-1 for i in range(self.n_case)]
            self.sz_viscs = [-1 for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_shear_zone(i)

    def import_directory(self, _dir, **kwargs):
        '''
        Import from a directory, look for groups and cases
        Inputs:
            _dir (str): directory to import
            kwargs
        '''
        assert(os.path.isdir(_dir))
        GroupP.CASE_SUMMARY.import_directory(self, _dir)

    def update_t660(self, i):
        '''
        update t660
        '''
        case_dir = self.ab_paths[i]
        # use the SLABPLOT class to read the slab_morph.txt file
        # and get the t660
        SlabPlot = SLABPLOT('foo')
        try:
            t660 = SlabPlot.GetTimeDepthTip(case_dir, 660e3)
        except SLABPLOT.SlabMorphFileNotExistError:
            t660 = -1.0
        self.t660s[i] = t660
        
    def update_shear_zone(self, i):
        '''
        Update shear zone properties
        ''' 
        case_dir = self.ab_paths[i]
        Visit_Options = self.VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        self.sz_viscs[i] = Visit_Options.options["SHEAR_ZONE_CONSTANT_VISCOSITY"]
        self.sz_depths[i] = Visit_Options.options["SHEAR_ZONE_CUTOFF_DEPTH"]
        self.sz_thicks[i] = Visit_Options.options["INITIAL_SHEAR_ZONE_THICKNESS"]

    def write_file(self, o_path):
        '''
        write file
        Inputs:
            o_path (str): path of output
        '''
        attrs_to_output = ['cases', 'steps', 'times', 'wallclocks', 't660s', "sz_thicks", "sz_depths"]
        headers = ['cases', 'steps', 'times (yr)', 'wallclocks (s)', 't660s (yr)', "sz_thicks (m)", "sz_depths (m)"]

        with open(o_path, 'w') as fout: 
            self.write(fout, attrs_to_output, headers) 
        
        print("%s: Write file %s" % (Utilities.func_name(), o_path))

    # todo_diagram
    def plot_diagram_t660(self, **kwargs):
        '''
        plot a diagram of t660
        '''
        fig_path = kwargs.get('fig_path', None)

        fig, ax = plt.subplots()
        x_name = "sz_thicks"
        y_name = "sz_depths"
        z_name = "t660s"

        # mesh data
        SZTT, SZDD, tt660 = self.mesh_data(x_name, y_name, z_name)
        
        # generate the plot
        ax.pcolormesh(SZTT, SZDD, tt660)

        # save figure
        if fig_path is not None:
            fig.savefig(fig_path)



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
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='A json file')
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
    elif (_commend in ['--json_option', '-jo']):
        # json options
        ShowJsonOption()
    elif _commend == 'create_group':
        CreateGroup(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()