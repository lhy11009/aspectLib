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

def test_reference_gamma():
    '''
    test for the reference case after iteration gamma
    '''
    source_case_dir = os.path.join(source_dir, "test_reference_gamma")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reference_gamma')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))
    wb_std_path = os.path.join(source_case_dir, 'case_std.wb')
    wb_path = os.path.join(output_dir, 'case.wb')
    assert(filecmp.cmp(wb_path, wb_std_path))


def test_fix_CDPT():
    '''
    test for
    fix the CDPT setup
    '''
    source_case_dir = os.path.join(source_dir, "test_fix_CDPT")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_CDPT')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_2d_consistent_reset_density():
    '''
    test for
    1. Adding reseting density utilites in a 2d consistent geometry
    '''
    source_case_dir = os.path.join(source_dir, "test_2d_consistent_reset_density")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_2d_consistent_reset_density')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_2d_consistent_particle():
    '''
    test for
    1. Usage of the particle method in a 2d consistent geometry
    '''
    source_case_dir = os.path.join(source_dir, "test_2d_consistent_particle")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_2d_consistent_particle')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_fix_reset_corner():
    '''
    test for
    1. coarsen the side of the plate in the minimum refinement function, 
    with a level assigned
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_fix_reset_corner")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_reset_corner')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_coarsen_side_level():
    '''
    test for
    1. coarsen the side of the plate in the minimum refinement function, 
    with a level assigned
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_coarsen_side_level")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_coarsen_side_level')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_fix_bd_temperature_2890():
    '''
    test for fixing boundary temperature with different values of box
    depth - a box depth of 2890e3, which is out of the bound of the da file
    assert:
    1. the prm rewrites the value of bottom temperature despite of a slighter larger range (2890e3)
    than the range give in the depth average file (~2850e3)
    '''
    source_case_dir = os.path.join(source_dir, "test_fix_bd_temperature_2890")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_bd_temperature_2890')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_fix_bd_temperature():
    '''
    test for fixing boundary temperature with different values of box
    depth
    assert:
    1. the prm file rewrite the value of bottom temperature 
    '''
    source_case_dir = os.path.join(source_dir, "test_fix_bd_temperature")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_fix_bd_temperature')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_coarsen_side():
    '''
    test for
    1. coarsen the side of the plate in the minimum refinement function
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_coarsen_side")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_coarsen_side')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_2d_consistent_sc():
    '''
    test for
    1. adjust geometry to 2890 km depth
    2. adjust length of the box by adding the trailing length of the sp plate
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_2d_consistent_sc")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_2d_consistent_sc')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

def test_adjust_geometry():
    '''
    test for
    1. adjust geometry to 2890 km depth
    2. adjust length of the box by adding the trailing length of the sp plate
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_adjust_geometry")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_adjust_geometry')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_reset_composition_viscosity():
    '''
    test for
    1. reset the composition viscosity
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_reset_composition_viscosity")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_reset_composition_viscosity')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_std.prm')
    prm_path = os.path.join(output_dir, 'case.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))


def test_adjust_solver_resolution():
    '''
    test for
    1. apply a new solver tolerance (iteration 20, tolerance = 0.001)
    2. adapting the particle method
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_adjust_solver_resolution")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_adjust_solver_resolution')
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


def test_schellart_particle():
    '''
    test for applying a vacancy at the side of the subducting plate
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_schellart_particle")
    json_path = os.path.join(source_case_dir, 'case0.json')
    output_dir = os.path.join(test_dir,'test_schellart_particle')
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


def test_composite_rheology():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_composite_rheology")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_composite_rheology')
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


def test_wb_new_ridge():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_wb_new_ridge")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_wb_new_ridge')
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


def test_side_plate():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_side_plate")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_side_plate')
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


def test_newtonian_setup():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_3d_newton")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_3d_newton')
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


def test_phase_transition():
    '''
    test for generating consistent inputs with the 2d cases
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_phase_transition")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_phase_transition')
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


def test_ov_transit():
    '''
    test for applying an transit section in the overiding plate
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_ov_transit")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_ov_transit')
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


def test_sp_side():
    '''
    test for applying a vacancy at the side of the subducting plate
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_sp_side")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_sp_side')
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


def test_sp_side():
    '''
    test for applying a vacancy at the side of the subducting plate
    Asserts:
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_sp_side_prescribed")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_sp_side_prescribed')
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


def test_fast_first_step_paraview():
    '''
    test for output the fast_first_step when paraview is used as the visualization software
    Asserts:
        the contants in the case_f.prm 
    '''
    # test 0
    source_case_dir = os.path.join(source_dir, "test_fast_first_step_paraview")
    json_path = os.path.join(source_case_dir, 'case.json')
    output_dir = os.path.join(test_dir,'test_fast_first_step_paraview')
    if os.path.isdir(output_dir):
        rmtree(output_dir)
    create_case_with_json(json_path, CASE, CASE_OPT)
    assert(os.path.isdir(output_dir))  # check case generation
    prm_std_path = os.path.join(source_case_dir, 'case_f_std.prm')
    prm_path = os.path.join(output_dir, 'case_f.prm')
    assert(filecmp.cmp(prm_path, prm_std_path))

# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

