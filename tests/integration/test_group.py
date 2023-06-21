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

def test_documentation_group_in_dir():
    '''
    test the implementation of documentation
    Asserts:
    '''
    group_dir = os.path.join(source_dir, "test_documentation_group_in_dir")
    assert(os.path.isdir(group_dir))
    GDoc = GDOC()
    GDoc.execute(group_dir, o_dir=test_local_dir)
    # assert something 
    mkd_file = os.path.join(test_local_dir, "documentation", "group_doc.mkd")
    mkd_file_std = os.path.join(group_dir, "group_doc_std.mkd")
    assert(os.path.isfile(mkd_file)) # assert file generation
    assert(filecmp.cmp(mkd_file, mkd_file_std)) # assert file contents


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

