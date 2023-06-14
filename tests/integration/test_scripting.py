# -*- coding: utf-8 -*-
r"""Test for foo.py

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
from shilofue.Scripting import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_scripting')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# make a new directory for test output
test_dir = os.path.join(".test", "test_scripting")
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)

####
# Utility functions
####

def test_explicit_import():
    '''
    explicitly import a function from a module
    Assert:
        the import contents agree with the standard one
    '''
    # todo_import
    # test 1
    ofile = os.path.join(test_dir, "test_explicit_import0.py")
    ofile_std = os.path.join(source_dir, "test_explicit_import_std0.py")
    assert(os.path.isfile(ofile_std)) # std file exists
    contents = ExplicitImport("tests.integration.fixtures.test_scripting.PlotStatistics", "STATISTICS_PLOT", "class")
    with open(ofile, 'w') as fout:
        fout.write(contents)
    # assert something 
    assert(filecmp.cmp(ofile, ofile_std))
    # test 2
    ofile1 = os.path.join(test_dir, "test_explicit_import1.py")
    ofile_std1 = os.path.join(source_dir, "test_explicit_import_std1.py")
    assert(os.path.isfile(ofile_std)) # std file exists
    contents = ExplicitImport("tests.integration.fixtures.test_scripting.PlotStatistics", "PlotFigure", "function")
    with open(ofile1, 'w') as fout:
        fout.write(contents)
    # assert something 
    assert(filecmp.cmp(ofile1, ofile_std1))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

