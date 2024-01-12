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
import pytest
import filecmp  # for compare file contents
import numpy as np
from shilofue.TwoDSubduction0.Group import *
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
test_local_dir = os.path.join(test_dir, 'test_TwoDSubduction_group')
if os.path.isdir(test_local_dir):
    rmtree(test_local_dir)
os.mkdir(test_local_dir)

source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_TwoDSubduction')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_case_summary():
    '''
    test case summary
    '''
    group_dir = os.path.join(source_dir, 'test_case_summary')
    
    # test 1: import directory
    Case_Summary = CASE_SUMMARY()
    o_path = os.path.join(test_local_dir, 'case_summary1.txt')
    o_path_std = os.path.join(group_dir, 'case_summary_std1.txt')
    Case_Summary.import_directory(group_dir)
    Case_Summary.write_file_if_update(o_path)
    assert(filecmp.cmp(o_path, o_path_std))

    # test 2: import directory and calculate t660
    Case_Summary = CASE_SUMMARY()
    o_path = os.path.join(test_local_dir, 'case_summary2.txt')
    o_path_std = os.path.join(group_dir, 'case_summary_std.txt')
    Case_Summary.import_directory(group_dir, actions=['t660'])
    Case_Summary.write_file_if_update(o_path)
    assert(filecmp.cmp(o_path, o_path_std))


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

