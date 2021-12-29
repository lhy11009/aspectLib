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
# import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.TwoDSubduction0.Group import *
from shilofue.TwoDSubduction0.Cases import CASE, CASE_OPT
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'TwoDSubduction', 'test_group')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_create_group():
    '''
    A test for the interfaces and functions to create a new group
    Asserts:
        case (directory) and files are generated
    '''
    # test 1
    json_path = os.path.join(source_dir, 'test.json')
    group_opt = GROUP_OPT()
    group_opt.read_json(json_path)
    feature_opt = group_opt.values[1][0]
    group = GROUP(CASE, CASE_OPT)
    group.read_json_base(group_opt.get_base_json_path())
    if os.path.isdir(group_opt.get_output_dir()):
        rmtree(group_opt.get_output_dir())  # remove old directory
    group.create_group(*group_opt.to_create_group())
    dir_stds = ['test_group_SA20.0', 'test_group_SA40.0', 'test_group_SA60.0']
    for dir_std in dir_stds:
        case_dir = os.path.join(group_opt.get_output_dir(), dir_std)
        print(case_dir)  # case directory
        assert(os.path.isdir(case_dir)) # case generation
        prm_path = os.path.join(case_dir, 'case.prm')
        assert(os.path.isfile(prm_path))
        wb_path = os.path.join(case_dir, 'case.wb')
        assert(os.path.isfile(wb_path))
    # test 2
    json_path = os.path.join(source_dir, 'test_EBA_CDPT.json')
    group_opt = GROUP_OPT()
    group_opt.read_json(json_path)
    feature_opt = group_opt.values[1][0]
    group = GROUP(CASE, CASE_OPT)
    group.read_json_base(group_opt.get_base_json_path())
    if os.path.isdir(group_opt.get_output_dir()):
        rmtree(group_opt.get_output_dir())  # remove old directory
    group.create_group(*group_opt.to_create_group())
    dir_stds = ['eba_cdpt_SA80.0_OA40.0', 'eba_cdpt_SA80.0_OA20.0', 'eba_cdpt_SA40.0_OA20.0']
    for dir_std in dir_stds:
        case_dir = os.path.join(group_opt.get_output_dir(), dir_std)
        print(case_dir)  # case directory
        assert(os.path.isdir(case_dir)) # case generation
        prm_path = os.path.join(case_dir, 'case.prm')
        assert(os.path.isfile(prm_path))
        wb_path = os.path.join(case_dir, 'case.wb')
        assert(os.path.isfile(wb_path))
    # check prm & wb file
    out_path = os.path.join(group_opt.get_output_dir(), 'eba_cdpt_SA80.0_OA40.0', 'case.prm')
    std_path = os.path.join(source_dir, 'eba_cdpt_SA80.0_OA40.0.prm')
    assert(filecmp.cmp(out_path, std_path))
    out_path = os.path.join(group_opt.get_output_dir(), 'eba_cdpt_SA80.0_OA40.0', 'case.wb')
    std_path = os.path.join(source_dir, 'eba_cdpt_SA80.0_OA40.0.wb')
    assert(filecmp.cmp(out_path, std_path))


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

