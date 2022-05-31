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
    First, substitute FOO with the name of the project
"""
import numpy as np
import sys, os, argparse
# import pathlib
import numpy as np
import shilofue.ParsePrm as ParsePrm
import shilofue.PlotVisit as PlotVisit
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, RunScripts

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
        Lib_FOO0_PlotVisit visit_options -i `pwd` -sr temperature.py\n\
            -sr: the relative path under the visit_script folder\n\
\n\
  - run script: \n\
        Lib_FOO0_PlotVisit run -i `pwd`\n\
\n\
  - run vtk scripts: \n\
        Lib_FOO0_PlotVisit vtk_options -i `pwd` -p FOO_MOW -s 25\n\
            -p operations, available operations are: \n\
                FOO_MOW: pull out tentative MOW from temperature: \n\
                FOO_SlabAnalysis: analyze slab morphology: \n\
        ")


class VISIT_OPTIONS(PlotVisit.VISIT_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
            last_step(list): plot the last few steps
        """
        # call function from parent
        PlotVisit.VISIT_OPTIONS.Interpret(self, **kwargs)
        idx = ParsePrm.FindWBFeatures(self.wb_dict, "Subducting plate")
        sub_plate_feature = self.wb_dict["features"][idx]
        sub_plate_extends = sub_plate_feature['coordinates'][1]
        self.options['TRENCH_EDGE_Y'] = sub_plate_extends[1] * 0.75


    def vtk_options(self, **kwargs):
        '''
        options of vtk scripts
        '''
        # call function from parent
        PlotVisit.VISIT_OPTIONS.vtk_options(self, **kwargs)


class PREPARE_RESULT_OPTIONS(PlotVisit.PREPARE_RESULT_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
        """
        # call function from parent
        PlotVisit.PREPARE_RESULT_OPTIONS.Interpret(self, **kwargs)


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
