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
from shutil import rmtree, copy2  # for remove directories

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
    print("base command:\n%s %s %s %s %s" % (executable, arg1, arg2, arg3, arg4))

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


def test_series_case_slurm():
    '''
    Here we test submitting a series of jobs to slurm.
    Each newer job would dependent on the completion of the previous job.
    '''
    executable = os.path.join(ASPECT_LAB_DIR, "bash_scripts", "parse_case.sh")

    file_sh =  os.path.join(big_source_dir, "eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr9_yd100", "job_p-billen.sh")
    file_prm =  os.path.join(big_source_dir, "eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr9_yd100", "case.prm")
    
    # assert 1: function of parse_slurm_file_series
    #   fileout generated and compared to the standard file
    target_dir = os.path.join(test_dir, "bash_submit_time_series")
    if os.path.isdir(target_dir):
        rmtree(target_dir)  # remove old outputs
    os.mkdir(target_dir)
    copy2(file_sh, target_dir)
    copy2(file_prm, target_dir)
    arg1 = "parse_slurm_file_series"
    arg2 = os.path.join(target_dir, "job_p-billen.sh")
    arg3 = "000001"
    arg4 = "0"
    completed_process = subprocess.run([executable, arg1, arg2, arg3, arg4], capture_output=False, text=True)
    print("base command:\n%s %s %s %s %s" % (executable, arg1, arg2, arg3, arg4))
    fileout = os.path.join(target_dir, "job_series_0.sh")
    assert(os.path.isfile(fileout))
    fileout_std = os.path.join(source_dir, "job_series_std.sh")
    assert(filecmp.cmp(fileout, fileout_std))

    # assert 2: function of submit_time_series
    arg1 = "submit_time_series"
    arg2 = target_dir
    arg3 = "job_p-billen.sh"
    arg4 = "3" # number of cases in the series
    arg5 = "1" # test_only, no case submitted to the system
    completed_process = subprocess.run([executable, arg1, arg2, arg3, arg4, arg5], capture_output=False, text=True)
    print("base command:\n%s %s %s %s %s %s" % (executable, arg1, arg2, arg3, arg4, arg5))
    fileout = os.path.join(target_dir, "job_series_1.sh")
    assert(os.path.isfile(fileout))
    fileout_std = os.path.join(source_dir, "job_series_1_std.sh")
    assert(filecmp.cmp(fileout, fileout_std))
    fileout = os.path.join(target_dir, "job_series_2.sh")
    assert(os.path.isfile(fileout))
    fileout = os.path.join(target_dir, "job_series_3.sh")
    assert(os.path.isfile(fileout))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

