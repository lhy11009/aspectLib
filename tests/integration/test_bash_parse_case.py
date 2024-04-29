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


from http.server import executable
import os
# import pytest
import filecmp  # for compare file contents
import subprocess

from shilofue.TwoDSubduction0.VtkPp import ASPECT_LAB_DIR
# import numpy as np
# import shilofue.Foo as Foo  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test/test_base_parse"
if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')
big_source_dir = os.path.join(ASPECT_LAB_DIR, "tests", "integration", "big_fixtures", "test_bash_parse")


def test_copy_case_output_by_step():
    '''
    test the option of copy_case_output_by_vtu_snapshot:
    give a case directory, copy the contents to a target directory and only include vtu files with a given vtu_snapshot number
    Asserts:
    '''
    executable = os.path.join(ASPECT_LAB_DIR, "bash_scripts", "parse_case.sh")

    arg1 = "copy_case_output_by_vtu_snapshot"
    arg2 =  os.path.join(big_source_dir, "eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr9_yd100") # this is the case directory to copy, so first check it's existence
    assert(os.path.isdir(arg2))
    arg3 = test_dir
    arg4 = "0"

    target_dir = os.path.join(test_dir, "eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr9_yd100")
    if os.path.isdir(target_dir):
        rmtree(target_dir)  # remove old outputs

    completed_process = subprocess.run([executable, arg1, arg2, arg3, arg4], capture_output=False, text=True)

    # assert something 
    assert(os.path.isdir(target_dir)) # assert new directory is generated
    prm_file = os.path.join(target_dir, "case.prm")  # assert files are copied to the new location
    assert(os.path.isfile(prm_file))
    wb_file = os.path.join(target_dir, "case.wb")
    assert(os.path.isfile(wb_file))
    sh_file = os.path.join(target_dir, "job_p-billen.sh")
    assert(os.path.isfile(sh_file))
    target_output_dir = os.path.join(target_dir, 'output') # assert the output directory is generated
    assert(os.path.isdir(target_output_dir)) # assert new directory is generated
    original_prm_file = os.path.join(target_output_dir, "original.prm")  # assert files are copied to the new location (output)
    assert(os.path.isfile(original_prm_file))
    statistic_file = os.path.join(target_output_dir, "statistics")
    assert(os.path.isfile(statistic_file))
    target_solution_dir=os.path.join(target_output_dir, "solution")
    assert(os.path.isdir(target_solution_dir)) # assert new directory is generated
    target_solution_pvtu_file=os.path.join(target_solution_dir, "solution-00000.pvtu") # assert files are copied to the new location (solution)
    assert(os.path.isfile(target_solution_pvtu_file))
    target_solution_visit_file=os.path.join(target_solution_dir, "solution-00000.visit")
    assert(os.path.isfile(target_solution_visit_file))
    target_solution_vtu_file=os.path.join(target_solution_dir, "solution-00000.0002.vtu")
    assert(os.path.isfile(target_solution_vtu_file))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

