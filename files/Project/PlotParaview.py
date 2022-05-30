# -*- coding: utf-8 -*-
r"""Scripts for working with paraview

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
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# locate import
from shilofue.FOO0.PlotVisit import VISIT_OPTIONS

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
  - prepare paraview scripts, this generates paraview scripts for a specific case: \n\
\n\
    python -m shilofue.PlotParaview paraview_options -i /home/lochy/ASPECT_PROJECT/ThDSubduction/s07_cases_1/s07_bl4000.0_bw4000.0_sw1000.0 -sr ThDSubduction/slab.py\n\
\n\
        ")

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
    parser.add_argument('-sr', '--script', type=str,
                        default='temperature.py',
                        help='script')
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
    elif _commend == 'paraview_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        case_dir = arg.inputs
        Paraview_Options = VISIT_OPTIONS(case_dir)
        # call function
        Paraview_Options.Interpret()
        # ofile = os.path.join('visit_scripts', 'slab_sph.py')
        ofile = os.path.join('paraview_scripts', os.path.basename(arg.script))
        paraview_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', arg.script)
        paraview_base_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', 'base.py')  # base.py : base file
        Paraview_Options.read_contents(paraview_base_script, paraview_script)  # this part combines two scripts
        Paraview_Options.substitute()  # substitute keys in these combined file with values determined by Interpret() function
        ofile_path = Paraview_Options.save(ofile, relative=True)  # save the altered script
        pass
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()