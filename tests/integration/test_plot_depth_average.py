# -*- coding: utf-8 -*-
r"""Tests for PlotDepthAverage.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest tests/integration/test_plot_depth_average.py 

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
import json
import pytest
# import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.Utilities import UNITCONVERT
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_depth_average')


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_plot_depth_average():
    '''
    A test on ploting depth average results
    Assertions:
        Errors are raised when inputs are not correct;
        The Inputs data are correct;
    '''
    _ofile_route = os.path.join(test_dir, 'DepthAverage.pdf')
    _ofile = os.path.join(test_dir, 'DepthAverage_t0.00000000e+00.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    # test_file = 'fixtures/statistics'
    test_file = os.path.join(source_dir, 'depth_average.txt')
    assert(os.access(test_file, os.R_OK))

    # use a json file
    json_file = os.path.join(source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('DepthAverage', {})

    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert, options=plot_options)
    # test error handling of key word time
    with pytest.raises(TypeError) as _excinfo:
        DepthAverage(test_file, fileout=_ofile_route, time='foo')
    assert(r'type of time' in str(_excinfo.value))
    # similar, but this time error is on an entry
    with pytest.raises(TypeError) as _excinfo:
        DepthAverage(test_file, fileout=_ofile_route, time=['foo'])
    assert(r'type of values in time' in str(_excinfo.value))
    # generate a successful output with time = 0.0
    DepthAverage(test_file, fileout=_ofile_route, time=0.0)
    assert(DepthAverage.time_step_length == 50)
    assert(DepthAverage.time_step_indexes[-1][-1] == 376)
    assert(abs(DepthAverage.time_step_times[0]-0.0) < 1e-6)
    assert(abs(DepthAverage.time_step_times[-1]-2.63571e+06)/2.63571e+06 < 1e-6)

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

