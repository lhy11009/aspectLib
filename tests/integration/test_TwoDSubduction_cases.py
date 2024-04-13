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


def test_sz_same_composition():
    '''
    values in the CDPT clapeyron slope
    '''
    source_case_dir = os.path.join(source_dir, "test_sz_same_composition")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_sz_same_composition')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_CDPT_clapeyron_slope():
    '''
    values in the CDPT clapeyron slope
    '''
    source_case_dir = os.path.join(source_dir, "test_CDPT_clapeyron_slope")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_CDPT_clapeyron_slope')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_reference_dry_olivine():
    '''
    The reference case after iteration gamma
    '''
    source_case_dir = os.path.join(source_dir, "test_reference_dry_olivine")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reference_dry_olivine')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_perple_X_table():
    '''
    The reference case after iteration gamma
    '''
    source_case_dir = os.path.join(source_dir, "test_perple_X_table")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_perple_X_table')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_reference_gamma_box():
    '''
    The reference case after iteration gamma
    '''
    source_case_dir = os.path.join(source_dir, "test_reference_gamma_box")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reference_gamma_box')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_reference_gamma():
    '''
    The reference case after iteration gamma
    '''
    source_case_dir = os.path.join(source_dir, "test_reference_gamma")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reference_gamma')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_fix_CDPT():
    '''
    fix the implementation of the CDPT phase transitions
    '''
    source_case_dir = os.path.join(source_dir, "test_fix_CDPT")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_CDPT')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))
    pass


def test_reference_case():
    '''
    to fix the conditions in a consistent model with the 3d models
    '''
    source_case_dir = os.path.join(source_dir, "test_reference_case")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reference_case')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_use_particle():
    '''
    to fix the conditions in a consistent model with the 3d models
    '''
    source_case_dir = os.path.join(source_dir, "test_use_particle")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_use_particle')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_remove_ov_comp():
    '''
    to fix the conditions in a consistent model with the 3d models
    '''
    source_case_dir = os.path.join(source_dir, "test_remove_ov_comp")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_remove_ov_comp')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))

def test_fix_3d_consistent():
    '''
    to fix the conditions in a consistent model with the 3d models
    '''
    source_case_dir = os.path.join(source_dir, "test_fix_3d_consistent")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_3d_consistent')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_box_dbg():
    '''
    to fix the issue in the wb file for box geometry
    '''
    source_case_dir = os.path.join(source_dir, "test_box_dbg")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_box_dbg')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_peierls_rheology_idrissi16():
    '''
    todo_idr
    test the idrissi 16 peierls rheology
    '''
    source_case_dir = os.path.join(source_dir, "peierls_rheology_idrissi16")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir, 'peierls_rheology_idrissi16')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_mantle_rheology_coh():
    '''
    apply a different Coh value in the mantle rheology
    '''
    source_case_dir = os.path.join(source_dir, "mantle_rheology_coh")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'mantle_rheology_coh')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))

def test_strong_slab_core():
    '''
    apply a strong core of the slab
    '''
    source_case_dir = os.path.join(source_dir, "strong_slab_core")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'strong_slab_core')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))

def test_sz_nl_visc_cpcl():
    '''
    Use the implementation for shear zone thickness
    also couple the viscosity to the eclogite transition
    '''
    source_case_dir = os.path.join(source_dir, "sz_nl_visc_cpcl")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_nl_visc_cpcl')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))

def test_trailing_edge_distance():
    '''
    Test with 2 layers in the crust
    '''
    source_case_dir = os.path.join(source_dir, "test_trailing_edge_distance")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_trailing_edge_distance')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))
    pass

def test_two_layer_crust():
    '''
    Test with 2 layers in the crust
    '''
    source_case_dir = os.path.join(source_dir, "test_2_layer")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_2_layer')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))
    pass


def test_fix_before_group1_old_sp():
    '''
    Tests before I fix the group of cases without the Peierls creep
    '''
    source_case_dir = os.path.join(source_dir, "fix_before_group1_old_sp")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'fix_before_group1_old_sp')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
#    assert(os.path.isdir(output_dir))  # check case generation
#    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
#    prm_path = os.path.join(output_dir, 'case.prm')
#    assert(filecmp.cmp(prm_path, prm_std_path))
#    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
#    wb_path = os.path.join(output_dir, 'case.wb')
#    assert(filecmp.cmp(wb_path, wb_std_path))


def test_bd_v():
    '''
    todo_bd
    test the boundary velcoty condition of "top prescrbed with bottom right open"
    '''
    source_case_dir = os.path.join(source_dir, "bd_v")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'bd_v')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


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

def test_mantle_rheology():
    '''
    Adjust the mantle rheology
    '''
    source_case_dir = os.path.join(source_dir, "mantle_rheology")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'mantle_rheology')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_shear_zone_strength():
    '''
    Use a stress dependent rheology in the shear zone
    '''
    source_case_dir = os.path.join(source_dir, "basalt_strengh_profile")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'basalt_strengh_profile')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_shear_zone_constant_viscosity():
    '''
    Use a constant viscosity in the shear zone
    '''
    source_case_dir = os.path.join(source_dir, "sz_constant_viscosity")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_constant_viscosity')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_wb_new_ridge_implementation():
    '''
    Use the new implemnetation for ridge coordinates in the world builder
    '''
    source_case_dir = os.path.join(source_dir, "new_ridge_implementation")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'new_ridge_implementation')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_wb_new_ridge_implementation_update():
    '''
    Use the new implemnetation for ridge coordinates in the world builder;
    also test the update file functionality of the create_case_with_json function.
    '''
    source_case_dir = os.path.join(source_dir, "new_ridge_implementation_update")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'new_ridge_implementation_update')
    output_dir_tmp = os.path.join(test_dir,'new_ridge_implementation_update_tmp')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    if os.path.isdir(output_dir_tmp):
        rmtree(output_dir_tmp)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))
    # now update on this existing case
    json_path = os.path.join(source_case_dir, 'case1.json')
    create_case_with_json(json_path, CASE, CASE_OPT, update=True, force_update=True)  # create case
    assert(os.path.isdir(os.path.join(output_dir, "update_00")))  # assert the update catalog is generated
    wb_std_path = os.path.join(source_case_dir, 'case_1_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(prm_path, prm_std_path)) # case.prm is unchanged
    assert(filecmp.cmp(wb_path, wb_std_path)) # case.wb is updated
    change_log_path = os.path.join(output_dir, "update_00", "change_log")
    assert(os.path.isfile(change_log_path))
    change_log_path_std = os.path.join(source_case_dir, "change_log")
    assert(filecmp.cmp(change_log_path, change_log_path_std)) # case.wb is updated
    # update again, with changes in the slurm file
    json_path = os.path.join(source_case_dir, 'case2.json')
    create_case_with_json(json_path, CASE, CASE_OPT, update=True, force_update=True)  # create case
    assert(os.path.isdir(os.path.join(output_dir, "update_01")))  # assert the update catalog is generated
    change_log_path = os.path.join(output_dir, "update_01", "change_log")
    assert(os.path.isfile(change_log_path))
    change_log_path_std = os.path.join(source_case_dir, "change_log_1")
    assert(filecmp.cmp(change_log_path, change_log_path_std)) # case.wb is updated


def test_sz_thickness():
    '''
    Use the implementation for shear zone thickness
    '''
    source_case_dir = os.path.join(source_dir, "sz_thickness")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_thickness')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_slurm_options():
    '''
    Use the implementation for slurm options
    '''
    source_case_dir = os.path.join(source_dir, "slurm_options")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'slurm_options')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    slurm_std_path = os.path.join(source_case_dir, 'job_p-billen_std.sh')
    slurm_path = os.path.join(output_dir, 'job_p-billen.sh')
    assert(filecmp.cmp(slurm_path, slurm_std_path))
    slurm_std_path_1 = os.path.join(source_case_dir, 'job_high2_std.sh')
    slurm_path_1 = os.path.join(output_dir, 'job_high2.sh')
    assert(filecmp.cmp(slurm_path_1, slurm_std_path_1))


def test_branch():
    '''
    Use the implementation for slurm options
    '''
    source_case_dir = os.path.join(source_dir, "test_branch")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_branch')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    slurm_std_path_1 = os.path.join(source_case_dir, 'job_high2_std.sh')
    slurm_path_1 = os.path.join(output_dir, 'job_high2.sh')
    assert(filecmp.cmp(slurm_path_1, slurm_std_path_1))


def test_old_sp_age():
    '''
    Use the implementation for shear zone thickness
    '''
    source_case_dir = os.path.join(source_dir, "old_sp_plate")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'old_sp_plate')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_sz_nl_visc():
    '''
    Use the implementation for shear zone thickness
    '''
    source_case_dir = os.path.join(source_dir, "sz_nl_visc")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_nl_visc')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_sz_ef():
    '''
    Use the embeded fault implementation for shear zone
    '''
    source_case_dir = os.path.join(source_dir, "sz_ef")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_ef')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))
    particle_std_path = os.path.join(source_case_dir, 'particle_std.dat')
    particle_path = os.path.join(output_dir, 'particle_file', 'particle.dat')
    assert(os.path.isfile(particle_path))
    assert(filecmp.cmp(particle_path, particle_std_path))


def test_sz_ef_feature_surface():
    '''
    Use the embeded fault implementation for shear zone, with the implementation of the 
    World Builder feature surface
    '''
    source_case_dir = os.path.join(source_dir, "sz_ef_feature_surface")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'sz_ef_feature_surface')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_eclogite_lowP():
    '''
    Test a setup for the eclogite transition that matches the mineral phase transtions
    '''
    source_case_dir = os.path.join(source_dir, "test_eclogite_lowP")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_eclogite_lowP')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)  # create case
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_0_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_0_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))



    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

