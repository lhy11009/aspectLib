# -*- coding: utf-8 -*-
r"""Test for foo.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_affinity.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
# import pytest
import filecmp  # for compare file contents
# import numpy as np
from shilofue.AffinityTest import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test/test_affinity"
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)

source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_affinity')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_affinity_base():
    '''
    (description)
    Asserts:
    '''
    test_source_dir = os.path.join(source_dir, "test_affinity_base")
    server = "peloton-rome"
    o_dir = os.path.join(test_dir, "test_affinity_base")
    os.mkdir(o_dir)
    tasks_per_node = 128
    refinement_levels = [2, 3, 4, 5]
    openmpi = "4.1.0"
    build_directory = "master_TwoD"
    prm_path = os.path.join(ASPECT_LAB_DIR, "files/AffinityTest/spherical_shell_expensive_solver.prm")
    assert(os.path.isfile(prm_path))
    slurm_base_path = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "220810", "job_p-billen.sh")
    # test 1: rene's affinity test 
    Affinity = AFFINITY(o_dir, prm_path, slurm_base_path, server, tasks_per_node, refinement_levels, -1, "block AMG", "", openmpi=openmpi, build_directory=build_directory)
    Affinity()
    # check the parent
    case_parent_dir = os.path.join(o_dir, "tmp", "peloton-rome-128tasks-socket-openmpi-4.1.0")
    assert(os.path.isdir(case_parent_dir))
    # check one child - 128 cpus, refinement level 4
    prm_128_4_path = os.path.join(case_parent_dir, "input_128_4_1", "case.prm")
    assert(os.path.isfile(prm_128_4_path))
    prm_128_4_path_std = os.path.join(test_source_dir, "input_128_4_1_std.prm")
    assert(filecmp.cmp(prm_128_4_path, prm_128_4_path_std))
    prm_128_4_slurm_path = os.path.join(case_parent_dir, "input_128_4_1", "job.sh")
    assert(os.path.isfile(prm_128_4_slurm_path))
    prm_128_4_slurm_path_std = os.path.join(test_source_dir, "input_128_4_1_std.sh")
    assert(filecmp.cmp(prm_128_4_slurm_path, prm_128_4_slurm_path_std))


def test_affinity_json():
    '''
    same tests as before, input from a json file
    '''
    test_source_dir = os.path.join(source_dir, "test_affinity_base")
    json_path = os.path.join(test_source_dir, "affinity_test.json")
    assert(os.path.isfile(json_path))
    create_tests_with_json(json_path, AFFINITY, AFFINITY_OPT)
    o_dir = os.path.join(test_dir, "test_affinity_base_1")
    # check the parent
    case_parent_dir = os.path.join(o_dir, "tmp", "peloton-rome-128tasks-socket-openmpi-4.1.0")
    assert(os.path.isdir(case_parent_dir))
    # check one child - 128 cpus, refinement level 4
    prm_128_4_path = os.path.join(case_parent_dir, "input_128_4_1", "case.prm")
    assert(os.path.isfile(prm_128_4_path))
    prm_128_4_path_std = os.path.join(test_source_dir, "input_128_4_1_std.prm")
    assert(filecmp.cmp(prm_128_4_path, prm_128_4_path_std))
    prm_128_4_slurm_path = os.path.join(case_parent_dir, "input_128_4_1", "job.sh")
    assert(os.path.isfile(prm_128_4_slurm_path))
    prm_128_4_slurm_path_std = os.path.join(test_source_dir, "input_128_4_1_std.sh")
    assert(filecmp.cmp(prm_128_4_slurm_path, prm_128_4_slurm_path_std))


def test_affinity_json_project():
    '''
    same tests as before, input from a json file
    This time from a project configureation
    '''
    test_source_dir = os.path.join(source_dir, "test_affinity_project")
    json_path = os.path.join(test_source_dir, "affinity_test.json")
    assert(os.path.isfile(json_path))
    create_tests_with_json(json_path, AFFINITY, AFFINITY_OPT)
    o_dir = os.path.join(test_dir, "test_affinity_project")
    # check the parent
    case_parent_dir = os.path.join(o_dir, "tmp", "stampede2-48tasks-socket-bGMG")
    assert(os.path.isdir(case_parent_dir))
    # check the prm file
    prm_file_1008_5_1 = os.path.join(case_parent_dir, "input_1008_5_1", "case.prm")
    prm_file_1008_5_1_std = os.path.join(test_source_dir, "case_1008_5_1_std.prm")
    assert(filecmp.cmp(prm_file_1008_5_1, prm_file_1008_5_1_std))


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))
