# -*- coding: utf-8 -*-
r"""Test for LatentHeatBK0/Cases.py

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
from shilofue.LatentHeatBK0.Cases import *
from shilofue.Cases import create_case_with_json
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test/LatentHeatBK_cases"
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', "test_LatentHeatBK")
source_dir1 = os.path.join(ASPECT_LAB_DIR, "files", "LatentHeatBK")

def test_layered():
    '''
    test the configure_wb function of class CASE
    here we use a different we to deploy tests than in the TwoDSubduction cases.
    I directly use the file saved under the "files" folder. But in order to do
    this, I need also fix the entries when calling the "create_case_with_json" function.
    In this function, I tested the cases with layered structures
    Asserts:
        1. prm file contains the right parameters
        2. cases are created successfully
    '''
    # test 1: test creating layered temperature structure
    source_case_dir = os.path.join(source_dir1, "03152022")
    source_case_dir_std = os.path.join(source_dir, "layered")
    json_path = os.path.join(source_case_dir, 'case_layered.json')
    output_dir = os.path.join(test_dir,'layered')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT, fix_case_name='layered',
    fix_base_dir=source_case_dir, fix_output_dir=test_dir)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_path=os.path.join(output_dir, 'case.prm')
    prm_std_path=os.path.join(source_case_dir_std, 'case_std.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_fix_velocity():
    '''
    In this function, I tested the cases where the vertical velocity is fixed
    Asserts:
        1. prm file contains the right parameters
        2. cases are created successfully
    '''
    source_case_dir = os.path.join(source_dir1, "03152022")
    source_case_dir_std = os.path.join(source_dir, "fix_velocity")
    json_path = os.path.join(source_case_dir_std, 'case_fix_velocity.json')
    output_dir = os.path.join(test_dir,'fix_velocity')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT, fix_case_name='fix_velocity',
    fix_base_dir=source_case_dir, fix_output_dir=test_dir)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_path=os.path.join(output_dir, 'case.prm')
    prm_std_path=os.path.join(source_case_dir_std, 'case_std.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))