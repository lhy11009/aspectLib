# -*- coding: utf-8 -*-
r"""Test for PostHefesto.py

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
import filecmp  # for compare file contents
import numpy as np
from shilofue.Plot import LINEARPLOT
from shilofue.PostHefesto import * # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_hefesto')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_read_n_output():
    '''
    test read in and output hefesto lookup table
    Asserts:
        output file is generated
    '''
    # input file
    filein = os.path.join(source_dir, 'lookup_table.txt')
    assert(os.path.isfile(filein))
    # output path
    fileout = os.path.join(test_dir, 'lookup_table_0.txt')
    if os.path.isfile(fileout):
        os.remove(fileout)
    # call processfunction
    Hefesto = HEFESTO()
    Hefesto.read_table(filein)
    Hefesto.ProcessHefesto(fileout)
    # assert something 
    assert(os.path.isfile(fileout))


def test_read_dimensions():
    '''
    test reading dimension information
    Asserts:
        read in the correct information
    '''
    # read data
    filein = os.path.join(source_dir, 'lookup_table.txt')
    Plotter = LINEARPLOT('hefesto', {})
    Plotter.ReadHeader(filein)
    Plotter.ReadData(filein)
    col_P = Plotter.header['Pi']['col']
    min, delta, number = ReadFirstDimension(Plotter.data[:, col_P])
    # check results
    tolerance = 1e-6
    assert(abs(min - 0.0) < tolerance)
    assert(abs(delta - 0.01) / 0.01 < tolerance)
    assert(abs(number - 4) / 4 < tolerance)
    # second dimension
    col_T = Plotter.header['Ti']['col']
    min, delta, number = ReadSecondDimension(Plotter.data[:, col_T])
    assert(abs(min - 800.0) / 800.0 < tolerance)
    assert(abs(delta - 1.0) / 1.0 < tolerance)
    assert(abs(number - 2) / 2 < tolerance)


        

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

