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
from shilofue.ThDSubduction0.Cases import *
from shilofue.Cases import create_case_with_json

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_schellart07():
    '''
    (description)
    Asserts:
    '''
    # test 0: change the size of the box
    source_case_dir = os.path.join(source_dir, "test_schellart_07")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_schellart_07')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))
    # test 1: change the size of the plate
    source_case_dir = os.path.join(source_dir, "test_schellart_07")
    json_path = os.path.join(source_case_dir, 'case1.json')
    output_dir = os.path.join(test_dir,'test_schellart_07_1')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case1_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case1_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))
    # test 2
    source_case_dir = os.path.join(source_dir, "test_schellart_07")
    json_path = os.path.join(source_case_dir, 'case2.json')
    output_dir = os.path.join(test_dir,'test_schellart_07_2')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case2_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_schellart07_newtonian():
    '''
    test for the newtonian cases created after the schellart model
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_schellart_07_newtonian")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_schellart_07_newtonian')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_2d_consistent():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_2d_consistent")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_2d_consistent')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

