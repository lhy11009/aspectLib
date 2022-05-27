# -*- coding: utf-8 -*-
r"""Test for TwoDSubduction0/Cases.py

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
from shilofue.TwoDSubduction0.Cases import *
from shilofue.Cases import create_case_with_json
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test/TwoDSubduction_cases"
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', "test_TwoDSubduction")

def test_wb_setup():
    '''
    test the configure_wb function of class CASE
    Asserts:
        1. wb file contains the right parameters
        2. prm file contains the right parameters
        3. cases are created successfully
    '''
    source_wb_dir = os.path.join(source_dir, "wb_setup")
    prm_file = os.path.join(source_wb_dir, 'case.prm')
    wb_file = os.path.join(source_wb_dir, 'case.wb')
    assert(os.access(prm_file, os.R_OK))
    assert(os.access(wb_file, os.R_OK))
    # test 0, chunk geometry and use transit ov plate
    json_file = os.path.join(source_wb_dir, 'configure.json')
    assert(os.access(json_file, os.R_OK))
    Case = CASE('wb_setup', prm_file, True, wb_inputs=wb_file)
    Case_Opt = CASE_OPT()
    Case_Opt.read_json(json_file)
    Case_Opt.check()
    Case.configure_prm(*Case_Opt.to_configure_prm())
    Case.configure_wb(*Case_Opt.to_configure_wb())
    # check the parameters
    wb_dict = Case.wb_dict
    assert(len(wb_dict['features']) == 5)  # this has 5 features
    i0 = ParsePrm.FindWBFeatures(wb_dict,'Overiding plate 1')  # transit plate
    assert(wb_dict['features'][i0]['coordinates'] == \
        [[35.972864236749224, -5.0], [35.972864236749224, 5.0],\
            [41.368793872261605, 5.0], [41.368793872261605, -5.0]]) # position
    # create new case
    Case.create(test_dir)
    case_dir = os.path.join(test_dir, 'wb_setup')
    case_prm_file = os.path.join(case_dir, 'case.prm')
    case_wb_file = os.path.join(case_dir, 'case.wb')
    assert(os.path.isfile(case_prm_file)) # assert files generated
    assert(os.path.isfile(case_wb_file))
    # test 1, chunk geometry and doesn't use transit ov plate
    json_file = os.path.join(source_wb_dir, 'configure_1.json')
    assert(os.access(json_file, os.R_OK))
    Case = CASE('wb_setup_1', prm_file, True, wb_inputs=wb_file)
    Case_Opt = CASE_OPT()
    Case_Opt.read_json(json_file)
    Case_Opt.check()
    Case.configure_prm(*Case_Opt.to_configure_prm())
    Case.configure_wb(*Case_Opt.to_configure_wb())
    # check the parameters
    wb_dict = Case.wb_dict
    assert(len(wb_dict['features']) == 4)  # this has 4 features
    prm_dict = Case.idict
    assert(prm_dict['Prescribed temperatures']['Temperature function']['Function constants'] == \
        "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=1.0647e+00,\\\n                             AGEOP=1.2614e+15, TS=2.730e+02, TM=1.6730e+03, K=1.000e-06, VSUB=1.5855e-09, PHILIM=1e-6")
    # create new case
    Case.create(test_dir)
    case_dir = os.path.join(test_dir, 'wb_setup_1')
    case_prm_file = os.path.join(case_dir, 'case.prm')
    case_wb_file = os.path.join(case_dir, 'case.wb')
    assert(os.path.isfile(case_prm_file)) # assert files generated
    assert(os.path.isfile(case_wb_file))
    # test 2, box geometry and doesn't use transit ov plate
    json_file = os.path.join(source_wb_dir, 'configure_2.json')
    assert(os.access(json_file, os.R_OK))
    Case = CASE('wb_setup_2', prm_file, True, wb_inputs=wb_file)
    Case_Opt = CASE_OPT()
    Case_Opt.read_json(json_file)
    Case_Opt.check()
    Case.configure_prm(*Case_Opt.to_configure_prm())
    Case.configure_wb(*Case_Opt.to_configure_wb())
    # check the parameters
    wb_dict = Case.wb_dict
    assert(len(wb_dict['features']) == 5)  # this has 5 features
    prm_dict = Case.idict
    # create new case
    Case.create(test_dir)
    case_dir = os.path.join(test_dir, 'wb_setup_2')
    case_prm_file = os.path.join(case_dir, 'case.prm')
    case_wb_file = os.path.join(case_dir, 'case.wb')
    assert(os.path.isfile(case_prm_file)) # assert files generated
    assert(os.path.isfile(case_wb_file))


def test_create_cases():
    '''
    test
    Asserts:
        1. wb file contains the right parameters
        2. prm file contains the right parameters
        3. cases are created successfully
    '''
    # test 1: test changing the ages of the plates
    source_case_dir = os.path.join(source_dir, "change_plate_ages")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'change_plate_ages_0')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    # test 2: test adjusting box width with plate age
    source_case_dir = os.path.join(source_dir, "change_plate_ages")
    json_path = os.path.join(source_case_dir, 'case1.json')
    output_dir = os.path.join(test_dir,'change_plate_ages_1')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_1_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_1_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_peierls_rheology():
    '''
    # test using the peierls rheology
    '''
    source_case_dir = os.path.join(source_dir, "peierls_rheology")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'peierls0')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case0_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case0_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_peierls_rheology_2_stages():
    '''
    # test using the peierls rheology
    '''
    source_case_dir = os.path.join(source_dir, "peierls_rheology_two_stage")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'peierls0_two_stage')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case0_std.prm')
    prm_std_path_1 = os.path.join(source_case_dir, 'case1_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case0_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    prm_path_1 = os.path.join(output_dir, 'case_1.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(prm_path_1, prm_std_path_1))
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_adjust_box():
    '''
    Adjust the width of the box
    '''
    source_case_dir = os.path.join(source_dir, "adjust_box_width")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'adjust_box0')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path))
    assert(filecmp.cmp(wb_path, wb_std_path))  


def test_3d_coarse_resolution():
    '''
    Adjust the width of the box
    '''
    # 3d_coarse
    source_case_dir = os.path.join(source_dir, "3d_coarse_resolution")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'3d_coarse_resolution0')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    prm_path = os.path.join(output_dir, 'case.prm')
    wb_path = os.path.join(output_dir, 'case.wb')
    # assert(filecmp.cmp(prm_path, prm_std_path))
    # assert(filecmp.cmp(wb_path, wb_std_path))  

    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

