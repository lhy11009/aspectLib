# -*- coding: utf-8 -*-
r"""Test for ParsePrm.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_parse_prm.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import shilofue.ParsePrm as ParsePrm
import os
# import pytest
import filecmp  # for compare file contents
# import numpy as np
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories


# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
test_dir = os.path.join(ASPECT_LAB_DIR, ".test")
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_prm')


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_parse_from_file():
    # test_file = 'fixtures/parse_test.prm'
    test_file = os.path.join(source_dir, 'case1.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = ParsePrm.ParseFromDealiiInput(fin)
    assert(inputs['Dimension'] == '2')
    assert(inputs['Use years in output instead of seconds'] == 'true')
    assert(inputs['End time'] == '40.0e6')
    assert(inputs['Additional shared libraries']
           == '/home/lochy/ASPECT_PROJECT/aspect_plugins/subduction_temperature2d/libsubduction_temperature2d.so, /home/lochy/ASPECT_PROJECT/aspect_plugins/prescribe_field/libprescribed_temperature.so')


def test_fast_first_step():
    '''
    (description)
    Asserts:
    '''
    _path = os.path.join(source_dir, 'case0.prm')
    std_path = os.path.join(source_dir, 'case0_std.prm')
    out_path = os.path.join(test_dir, 'fast_first_step.prm')
    inputs = ParsePrm.ReadPrmFile(_path)
    ParsePrm.FastZeroStep(inputs)
    ParsePrm.WritePrmFile(out_path, inputs)
    # assert something 
    assert(filecmp.cmp(out_path, std_path))


def test_UpperMantleRheologyViscoPlastic():
    '''
    test UpperMantleRheologyViscoPlastic
    Asserts:
        The two read-in dictionaries are identical to standard
    '''
    _path = os.path.join(source_dir, 'case2.prm')
    inputs = ParsePrm.ReadPrmFile(_path)
    diffusion_creep, dislocation_creep = ParsePrm.UpperMantleRheologyViscoPlastic(inputs)
    a = str(diffusion_creep)
    b = str(dislocation_creep)
    assert(a == "{'A': 8.1787e-17, 'd': 0.01, 'n': 1.0, 'm': 3.0, 'E': 300000.0, 'V': 6.9e-06}")
    assert(b == "{'A': 5.9076e-16, 'd': 0.01, 'n': 3.5, 'm': 0.0, 'E': 510000.0, 'V': 1.74e-05}")

def test_ParseFromSlurmBatchFile():
    '''
    test function ParseFromSlurmBatchFile
    Asserts:
        options from reading the file, respectively
    '''
    source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_prm')
    _path = os.path.join(source_dir, 'job_p-billen.sh')
    assert(os.path.isfile(_path))
    with open(_path, 'r') as fin:
        i_dict = ParsePrm.ParseFromSlurmBatchFile(fin)
    # assertion
    # 1. header
    assert(len(i_dict["header"]) == 1 and i_dict["header"][0] == "#!/bin/bash -l")
    # 2. configuration
    assert(i_dict["config"]["--threads-per-core"] == '2')
    # 3. load
    assert(len(i_dict["load"]) == 2 and 'openmpi/4.1.0-mpi-io' in i_dict["load"])
    # 4. command
    assert(len(i_dict["command"]) == 3 and i_dict["command"][2] == "case.prm")
    pass


def test_ParseToSlurmBatchFile():
    '''
    test function ParseToSlurmBatchFile
    '''
    # test 1 read and write
    source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_prm')
    _path = os.path.join(source_dir, 'job_p-billen.sh')
    o_path = os.path.join(test_dir, "slurm.sh")
    if os.path.isfile(o_path):
        os.remove(o_path)
    o_std_path = os.path.join(source_dir, "slurm_std.sh")
    with open(_path, 'r') as fin:
        i_dict = ParsePrm.ParseFromSlurmBatchFile(fin)
    with open(o_path, 'w') as fout:
        ParsePrm.ParseToSlurmBatchFile(fout, i_dict)
    assert(os.path.isfile(o_path))  # compare outputs
    assert(filecmp.cmp(o_path, o_std_path))
    # test 2: use class, read, change affinity and write
    source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_prm')
    _path = os.path.join(source_dir, 'job_p-billen.sh')
    o_path = os.path.join(test_dir, "slurm_test2.sh")
    if os.path.isfile(o_path):
        os.remove(o_path)
    o_std_path = os.path.join(source_dir, "slurm_std_test2.sh")
    SlurmOperator = ParsePrm.SLURM_OPERATOR(_path)
    SlurmOperator.SetAffinity(2, 64, 1, partition="high2")
    SlurmOperator(o_path)
    assert(os.path.isfile(o_path))  # compare outputs
    assert(filecmp.cmp(o_path, o_std_path))

def test_ParseToSlurmBatchFile_mpirun():
    '''
    test function ParseToSlurmBatchFile, and use the mpirun command
    for execution
    '''
    source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_prm')
    _path = os.path.join(source_dir, 'job_p-billen.sh')
    o_path = os.path.join(test_dir, "slurm_test3.sh")
    if os.path.isfile(o_path):
        os.remove(o_path)
    o_std_path = os.path.join(source_dir, "slurm_std_test3.sh")
    SlurmOperator = ParsePrm.SLURM_OPERATOR(_path)
    SlurmOperator.SetAffinity(2, 64, 1, partition="high2", use_mpirun=1, bind_to="socket")
    SlurmOperator(o_path)
    assert(os.path.isfile(o_path))  # compare outputs
    assert(filecmp.cmp(o_path, o_std_path))
