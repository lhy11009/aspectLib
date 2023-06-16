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
# import filecmp  # for compare file contents
import numpy as np
from shilofue.Scripting import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_parse_import_syntax():
    '''
    Asserts: module and object are parsed correctly
    '''
    # test 1
    line = "from shilofue.Plot import LINEARPLOT"
    module, _object = ParseImportSyntax(line)
    assert(module == "shilofue.Plot")
    assert(_object == ["LINEARPLOT"])
    # test 2: test the as syntax is eliminated
    line = "from shilofue.Plot import LINEARPLOT as foo"
    module, _object = ParseImportSyntax(line)
    assert(module == "shilofue.Plot")
    assert(_object == ["LINEARPLOT"])
    # test 3: assert import multiple objects
    line = "from shilofue.Plot import LINEARPLOT, foo"
    module, _object = ParseImportSyntax(line)
    assert(module == "shilofue.Plot")
    assert(_object == ["LINEARPLOT", "foo"])

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

