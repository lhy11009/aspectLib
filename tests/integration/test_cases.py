# -*- coding: utf-8 -*-
r"""Test for Cases.py

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
from shutil import rmtree  # for remove directoriesa
from shilofue.Cases import *

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def ConfigureFoo(Inputs, _config):
    """
    an example of configuation
    """
    return Inputs


def ConfigureFoo1(Inputs, _config):
    """
    another example of configuation, second option for renaming
    """
    return Inputs, "_foo"


def test_create_case():
    '''
    Use the interface defined in Cases.py. Take a inputs file, do a little multilation and create a new case
    Asserts:
        cases in created(files, contents, etc)
        assert prm file is generated
        assert prm file for fast running 0th step is generated
    '''
    prm_path = os.path.join(source_dir, 'case.prm')
    extra_path = os.path.join(source_dir, 'particle.dat')
    test_case = CASE('foo', prm_path, False)
    case_output_dir = os.path.join(test_dir, 'foo')
    if os.path.isdir(case_output_dir):
        rmtree(case_output_dir)
    test_case.configure(ConfigureFoo, 1)  # do nothing, test interface
    test_case.add_extra_file(extra_path)  # add an extra file
    # todo_intial
    test_case.create(test_dir, fast_first_step=1, test_initial_steps=(3, 1e4))
    # assert prm file is generated
    prm_output_path = os.path.join(case_output_dir, 'case.prm')
    prm_std_path = os.path.join(source_dir, 'case_std.prm')
    assert(os.path.isfile(prm_output_path))  # assert generated
    assert(filecmp.cmp(prm_output_path, prm_std_path))  # assert contents
    # assert prm file for fast running 0th step is generated
    prm_output_path = os.path.join(case_output_dir, 'case_f.prm')
    prm_std_path = os.path.join(source_dir, 'case_f_std.prm')
    assert(os.path.isfile(prm_output_path))  # assert generated
    assert(filecmp.cmp(prm_output_path, prm_std_path))  # assert contents
    # assert prm file for testing hte initial steps are generated
    prm_output_path = os.path.join(case_output_dir, 'case_ini.prm')
    prm_std_path = os.path.join(source_dir, 'case_ini_std.prm')
    assert(os.path.isfile(prm_output_path))  # assert generated
    assert(filecmp.cmp(prm_output_path, prm_std_path))  # assert contents
    # assert extra file is generated
    extra_output_path = os.path.join(case_output_dir, 'particle.dat')
    assert(os.path.isfile(extra_output_path))  # assert generated

    # test 2: renaming
    test_case.configure(ConfigureFoo1, 1, rename=True)  # do nothing, test interface
    test_case.create(test_dir)
    # renaming: add string '_foo'
    case_output_dir = os.path.join(test_dir, 'foo_foo')
    # assert prm file is generated
    prm_output_path = os.path.join(case_output_dir, 'case.prm')
    prm_std_path = os.path.join(source_dir, 'case_std.prm')
    assert(os.path.isfile(prm_output_path))  # assert generated
    assert(filecmp.cmp(prm_output_path, prm_std_path))  # assert contents
    # assert extra file is generated
    extra_output_path = os.path.join(case_output_dir, 'particle.dat')
    assert(os.path.isfile(extra_output_path))  # assert generated

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

