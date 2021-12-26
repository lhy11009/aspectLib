# -*- coding: utf-8 -*-
r"""to hold the CASE_OPTIONS class

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
from shilofue.PlotStatistics import STATISTICS_PLOT
from shilofue.ParsePrm import ParseFromDealiiInput

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

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


class CASE_OPTIONS(Utilities.CODESUB):
    """
    parse .prm file to a option file that bash can easily read
    This inherit from Utilities.CODESUB
    Attributes:
        _case_dir(str): path of this case
        _output_dir(str): path of the output
        _visit_file(str): path of the visit file
        options(dict): dictionary of key and value to output
        i_dict(dict): dictionary for prm file
        wb_dict(dict): dictionary for wb file
    """
    def __init__(self, case_dir):
        """
        Initiation
        Args:
            case_dir(str): directory of case
        """
        Utilities.CODESUB.__init__(self)
        # check directory
        self._case_dir = case_dir
        Utilities.my_assert(os.path.isdir(self._case_dir), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case directory - %s doesn\'t exist' % self._case_dir)
        self._output_dir = os.path.join(case_dir, 'output')
        Utilities.my_assert(os.path.isdir(self._output_dir), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case output directory - %s doesn\'t exist' % self._output_dir)
        self._visit_file = os.path.join(self._output_dir, 'solution.visit')
        Utilities.my_assert(os.access(self._visit_file, os.R_OK), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case visit file - %s cannot be read' % self._visit_file)
        # output dir
        self._output_dir = os.path.join(case_dir, 'output')
        if not os.path.isdir(self._output_dir):
            os.mkdir(self._output_dir)
        # img dir
        self._img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(self._img_dir):
            os.mkdir(self._img_dir)

        # get inputs from .prm file
        prm_file = os.path.join(self._case_dir, 'case.prm')
        Utilities.my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case prm file - %s cannot be read' % prm_file)
        with open(prm_file, 'r') as fin:
            self.idict = ParseFromDealiiInput(fin)
        
        # wb inputs:
        #   if there is a .wb file, it is loaded. Otherwise, just start with
        # a vacant dictionary.
        self.wb_dict = {}
        wb_file = os.path.join(self._case_dir, 'case.wb')
        if os.access(wb_file, os.R_OK):
            with open(wb_file, 'r') as fin:
                self.wb_dict = json.load(fin)

        # initiate a dictionary
        self.options = {}

        # initiate a statistic data
        self.Statistics = STATISTICS_PLOT('Statistics')
        self.statistic_file = os.path.join(self._output_dir, 'statistics')
        self.Statistics.ReadHeader(self.statistic_file)
        self.Statistics.ReadData(self.statistic_file)

        # horiz_avg
        self.horiz_avg_file = os.path.join(self._output_dir, "depth_average.txt")


    def Interpret(self):
        """
        Interpret the inputs, to be reloaded in children
        """
        pass
    
    def get_geometry(self):
        '''
        get the name of geomery
        '''
        return self.idict['Geometry model']['Model name']

    
    def save(self, _path, **kwargs):
        '''
        save contents to a new file
        Args:
            kwargs(dict):
                relative: use relative path
        '''
        use_relative_path = kwargs.get('relative', False)
        if use_relative_path:
            _path = os.path.join(self._case_dir, _path)
        o_path = Utilities.CODESUB.save(self, _path)
        return o_path

    def __call__(self, ofile, kwargs):
        """
        Call function
        Args:
            ofile(str): path of output
        """
        # interpret
        self.Interpret(kwargs)

        # open ofile for output
        # write outputs by keys and values
        with open(ofile, 'w') as fout:
            for key, value in self.options.items():
                fout.write("%s       %s\n" % (key, value))
        pass

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
    elif _commend == 'foo':
        # example:
        SomeFunction('foo')
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()