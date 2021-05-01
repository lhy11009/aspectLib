# -*- coding: utf-8 -*-
r"""Test for PlotRunTime.py

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
# import pytest
# import filecmp  # for compare file contents
# import numpy as np
import shilofue.PlotRunTime as PlotRunTime  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_run_time')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_plot_run_time():
    '''
    test read from a log file and then generate run time plot
    Asserts:
    '''
    ifile = os.path.join(source_dir, "log.txt")
    o_path = os.path.join(test_dir, "RunTime.png")
    if os.path.isfile(o_path):
        os.remove(o_path)
    PlotRunTime.PlotFigure(ifile, o_path) 

    # assert something 
    assert(os.path.isfile(o_path))  # assert figure is generated


def test_plot_run_time_newtion():
    '''
    test read from a log file and then generate run time plot of newton solver
    Asserts:
    '''
    ifile = os.path.join(source_dir, "log.txt")
    o_path_base = os.path.join(test_dir, "RunTimeNewton.png")
    o_path_std = os.path.join(test_dir, "RunTimeNewton_s000000.png")
    if os.path.isfile(o_path_std):
        os.remove(o_path_std)
    o_path = PlotRunTime.PlotNewtonSolver(ifile, o_path_base) 

    # assert something 
    assert(o_path == o_path_std)
    assert(os.path.isfile(o_path))  # assert figure is generated

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

