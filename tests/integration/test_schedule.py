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
from shilofue.Schedule import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'schedule')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_scheduler():
    '''
    test the implementation of the SCHEDULER class
    Asserts:
    '''
    o_dir = os.path.join(test_dir, 'test_scheduler')
    if not os.path.isdir(o_dir):
        os.mkdir(o_dir)
    # test 1
    o_local_file = os.path.join(o_dir, 'local.sh')
    o_remote_file = os.path.join(o_dir, 'remote.sh')
    o_local_file_std = os.path.join(source_dir, 'local_std.sh')
    o_remote_file_std = os.path.join(source_dir, 'remote_std.sh')
    if os.path.isfile(o_local_file):
        os.remove(o_local_file)  # remove old results
    if os.path.isfile(o_remote_file):
        os.remove(o_remote_file)  # remove old results
    assert(os.path.isfile(o_local_file_std))
    assert(os.path.isfile(o_remote_file_std))
    Scheduler = SCHEDULER()  # initiation
    new_task_local = ["foo", "foo1"]
    Scheduler.add_task(new_task_local)
    Scheduler.add_task("foo1", remote=True)
    Scheduler.add_task("foo2", remote=True)
    Scheduler(o_dir)
    # assert something 
    assert(filecmp.cmp(o_local_file, o_local_file_std))
    assert(filecmp.cmp(o_remote_file, o_remote_file_std))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

