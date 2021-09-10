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
from shilofue.PostHefesto import *
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'post_hefesto')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_process_hefesto_table():
    '''
    Test processing hefesto table
    Asserts:
    '''
    input_file = os.path.join(source_dir, "example_table")
    assert(os.path.isfile(input_file))  # assert there is an existing Hefesto table
    output_file = os.path.join(test_dir, "example_table_output")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "example_table_output_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    ProcessHefesto(input_file, output_file, 1, 1)  # process Hefesto output
    assert(os.path.isfile(output_file))  # assert there is output.
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))

    # assert something
    assert(True)


# notes

# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

