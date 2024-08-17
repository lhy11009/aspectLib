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
import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.Group import  *
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
test_local_dir = os.path.join(test_dir, "test_group")
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_group')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# make new test directory
if os.path.isdir(test_local_dir):
    rmtree(test_local_dir)
os.mkdir(test_local_dir)

def test_case_summary_latex():
    '''
    test writing a latex table
    '''
    # initiate
    Case_Summary = CASE_SUMMARY()

    # import a group 
    group_dir = os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1")
    assert(os.path.isdir(group_dir))
    Case_Summary.import_directory(group_dir)

    # write outputs
    o_file = os.path.join(test_local_dir, "case_summary.tex")
    o_file_std =  os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1", "ofile_std.tex")
    Case_Summary.write_file(o_file)
    assert(os.path.isfile(o_file))
    assert(os.path.isfile(o_file_std))
    assert(filecmp.cmp(o_file, o_file_std))


def test_case_summary():
    '''
    Test the CASE_SUMMARY class
    '''
    # initiate
    Case_Summary = CASE_SUMMARY()

    # import a group 
    group_dir = os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1")
    assert(os.path.isdir(group_dir))
    Case_Summary.import_directory(group_dir)

    # write outputs
    o_file = os.path.join(test_local_dir, "case_summary.txt")
    o_file_std =  os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1", "ofile_std")
    Case_Summary.write_file(o_file)
    assert(os.path.isfile(o_file))
    assert(os.path.isfile(o_file_std))
    assert(filecmp.cmp(o_file, o_file_std))

    # read file
    # file is generated in the last test
    Case_Summary1 = CASE_SUMMARY()
    Case_Summary1.import_txt(o_file)
    assert(Case_Summary1.cases == ["eba3d_SA80.0_OA40.0_width61_GR4_AR4", "eba3d_SA80.0_OA40.0_width61_GR3_AR3_sc_1e23"])
    assert(Case_Summary1.wallclocks == ["221.0", "85100.0"])


def test_case_summary_csv():
    '''
    Test the CASE_SUMMARY class, output csv format
    '''
    # initiate
    Case_Summary = CASE_SUMMARY()

    # import a group 
    group_dir = os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1")
    assert(os.path.isdir(group_dir))
    Case_Summary.import_directory(group_dir)

    # write outputs
    o_file = os.path.join(test_local_dir, "case_summary.csv")
    o_file_std =  os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1", "ofile_std.csv")
    Case_Summary.write_file(o_file)
    assert(os.path.isfile(o_file))
    assert(os.path.isfile(o_file_std))
    assert(filecmp.cmp(o_file, o_file_std))

    # read file
    # file is generated in the last test
    Case_Summary1 = CASE_SUMMARY()
    Case_Summary1.import_file(o_file)
    assert(Case_Summary1.cases == ["eba3d_SA80.0_OA40.0_width61_GR4_AR4", "eba3d_SA80.0_OA40.0_width61_GR3_AR3_sc_1e23"])
    assert(Case_Summary1.wallclocks == ["221.0", "85100.0"])


def test_documentation_group_in_dir():
    '''
    test the implementation of documentation
    Asserts:
    '''
    project_dir = os.path.join(source_dir, "test_documentation_group_in_dir")
    assert(os.path.isdir(project_dir))
    GDoc = GDOC()
    GDoc.execute(project_dir, o_dir=test_local_dir)
    # assert something 
    mkd_file = os.path.join(test_local_dir, "documentation", "group_doc.mkd")
    mkd_file_std = os.path.join(project_dir, "group_doc_std.mkd")
    assert(os.path.isfile(mkd_file)) # assert file generation
    assert(filecmp.cmp(mkd_file, mkd_file_std)) # assert file contents
    latex_file = os.path.join(test_local_dir, "documentation", "group_doc.tex")
    latex_file_std = os.path.join(project_dir, "group_doc_std.tex")
    assert(os.path.isfile(latex_file)) # assert file generation
    assert(filecmp.cmp(latex_file, latex_file_std)) # assert file contents


def test_ReadBasicInfoGroup():
    '''
    test the function ReadBasicInfoGroup
    assert outputs
    '''
    group_dir = os.path.join(source_dir, "test_documentation_group_in_dir", "EBA_2d_consistent_1")
    assert(os.path.isdir(group_dir))
    case_list, step_list, time_list, wallclock_list = ReadBasicInfoGroup(group_dir)
    # assert outputs
    assert(case_list==['eba3d_SA80.0_OA40.0_width61_GR4_AR4', 'eba3d_SA80.0_OA40.0_width61_GR3_AR3_sc_1e23'])
    assert(step_list == [1, 130])
    assert(time_list ==  [100000.0, 3724260.0])
    assert(wallclock_list == [221.0, 85100.0])


# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

