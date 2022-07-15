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
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories
from shilofue.PhaseDiagram0.Cases import *
from shilofue.Cases import create_case_with_json

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'PhaseDiagram', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_CDPT():
    '''
    (description)
    Asserts:
    '''
    # test 0: change compositions in the CDPT model
    source_case_dir = os.path.join(source_dir, "test_CDPT")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_CDPT')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_bd_lsolver():
    '''
    test for cases with bd and linear solver
    '''
    # todo_bc
    # test 1: tangential bc
    source_case_dir = os.path.join(source_dir, "test_bd_lsolver")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_bd_lsolver')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    # test 2: tangential bc on top and bottom, no slip on both sides
    source_case_dir = os.path.join(source_dir, "test_bd_lsolver_ns")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_bd_lsolver_ns')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    # test 3: tangential bc on top and bottom, no slip on both sides; box geometry
    source_case_dir = os.path.join(source_dir, "test_bd_lsolver_ns_box")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_bd_lsolver_ns_box')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_bd_lsolver_3D():
    '''
    test for cases with bd and linear solver
    '''
    # test 1: tangential bc, for 3d
    source_case_dir = os.path.join(source_dir, "test_bd_lsolver_box_3D")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_bd_lsolver_box_3D')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    # test 2: tangential bc on top and bottom, no slip on both sides, for 3d
    source_case_dir = os.path.join(source_dir, "test_bd_lsolver_ns_box_3D")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_PhaseDiagram_bd_lsolver_ns_box_3D')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

