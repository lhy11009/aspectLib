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
# import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.TwoDSubduction0.Cases import *
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test/TwoDSubduction_cases"
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', "test_TwoDSubduction")

def test_wb_setup():
    '''
    (description)
    Asserts:
    '''
    source_wb_dir = os.path.join(source_dir, "wb_setup")
    prm_file = os.path.join(source_wb_dir, 'case.prm')
    wb_file = os.path.join(source_wb_dir, 'case.wb')
    json_file = os.path.join(source_wb_dir, 'configure.json')
    assert(os.access(prm_file, os.R_OK))
    assert(os.access(wb_file, os.R_OK))
    assert(os.access(json_file, os.R_OK))
    Case = CASE('wb_setup', prm_file, wb_inputs=wb_file)
    Case_Opt = CASE_OPT()
    Case_Opt.read_json(json_file)
    Case.configure_wb(*Case_Opt.to_configure_wb())
    Case.create(test_dir)
    case_dir = os.path.join(test_dir, 'wb_setup')
    case_prm_file = os.path.join(case_dir, 'case.prm')
    case_wb_file = os.path.join(case_dir, 'case.wb')
    assert(os.path.isfile(case_prm_file)) # assert files generated
    assert(os.path.isfile(case_wb_file))
    # assert the context of the generated file
    # assert something 
    assert(True)

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

