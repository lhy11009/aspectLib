# -*- coding: utf-8 -*-
r"""Test for PlotStatistics.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_foo.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
import json
# import pytest
# import filecmp  # for compare file contents
# import numpy as np
from shilofue.PlotStatistics import STATISTICS_PLOT  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_statistics')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_foo():
    '''
    (description)
    Asserts:
    '''

    # assert something 
    assert(True)


def test_plot_statistics():
    '''
    A test on ploting statistics results
    '''
    _ofile = os.path.join(test_dir, 'Statistics.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    test_file = os.path.join(source_dir, 'statistics')
    assert(os.access(test_file, os.R_OK))

    # use a json file
    json_file = os.path.join(source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('Statistics', {})

    # plot statistics ouput #####
    Statistics = STATISTICS_PLOT('Statistics', options=plot_options)
    Statistics(test_file, fileout=_ofile)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished


    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

