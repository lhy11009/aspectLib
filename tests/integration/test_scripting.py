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
import re
# import pytest
import filecmp  # for compare file contents
import numpy as np
from shilofue.Scripting import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_scripting')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# make a new directory for test output
test_dir = os.path.join(".test", "test_scripting")
if os.path.isdir(test_dir):
    rmtree(test_dir)
os.mkdir(test_dir)


def test_Utilities():
    '''
    test function read header
    '''
    file_path = os.path.join(source_dir, "AffinityTest.py")
    # file_path1 = os.path.join(source_dir, "Rheology.py")
    assert(os.path.isfile(file_path))
    # assert(os.path.isfile(file_path1))
    # test 1: ReadInHeader1 function
    headers = ReadInHeaders1(file_path)
    assert('import shilofue.Cases as CasesP' in headers)
    assert('import shilofue.TwoDSubduction0.Cases as CasesTwoDSubduction' in headers) 
    assert('import shilofue.ThDSubduction0.Cases as CasesThDSubduction' in headers)
    assert('import shilofue.ParsePrm as ParsePrm' in headers)
    # test 2: find imported modules from file contents
    with open(file_path, 'r') as fin:
        slines = fin.readlines()
    module, alias = ParseModuleObject("import shilofue.Cases as CasesP") # parse the module and alias
    objects = FindImportModule(module, alias, slines)
    assert(module=="shilofue.Cases")
    assert(alias=="CasesP")
    assert("create_case_with_json" in objects)
    # test 3: ExplicitImport function
    # assert: the contents to import
    module = "shilofue.Cases"
    _object = "create_case_with_json"
    alias = "CasesP"
    explicit_import_contents = ExplicitImport(module, _object, alias=alias)
    assert(re.match(".*CasesP_create_case_with_json.*", explicit_import_contents.split('\n')[0]))
    # todo_script
    # test 4: find imported modules from file contents recursively
#    with open(file_path1, 'r') as fin:
#        slines = fin.readlines()
#    module, alias, objects = FindImportModuleRecursive("import shilofue.PlotDepthAverage as PDAver", slines)
#    print("module: ", module)
#    print("alias: ", alias)
#    print("objects: ", objects)
    # assert(module=="shilofue.Cases")
    # assert(alias=="CasesP")
    # assert("create_case_with_json" in objects)


def test_explicit_import():
    '''
    explicitly import a function from a module
    Assert:
        the import contents agree with the standard one
    '''
    # test 1
    ofile = os.path.join(test_dir, "test_explicit_import0.py")
    ofile_std = os.path.join(source_dir, "test_explicit_import_std0.py")
    assert(os.path.isfile(ofile_std)) # std file exists
    contents = ExplicitImport("tests.integration.fixtures.test_scripting.PlotStatistics", "STATISTICS_PLOT", "class")
    with open(ofile, 'w') as fout:
        fout.write(contents)
    # assert something 
    assert(filecmp.cmp(ofile, ofile_std))
    # test 2
#    ofile1 = os.path.join(test_dir, "test_explicit_import1.py")
#    ofile_std1 = os.path.join(source_dir, "test_explicit_import_std1.py")
#    assert(os.path.isfile(ofile_std)) # std file exists
#    contents = ExplicitImport("tests.integration.fixtures.test_scripting.PlotStatistics", "PlotFigure", "function")
#    with open(ofile1, 'w') as fout:
#        fout.write(contents)
#    # assert something 
#    assert(filecmp.cmp(ofile1, ofile_std1))


def test_parse_header():
    '''
    Parse the header of the input file, figuring out
    what to import
    assert:
        the right modules and objects are parsed from a file
    '''
    ifile = os.path.join(source_dir, "PlotStatistics.py")
    assert(os.path.isfile(ifile))
    module_list, object_list = ParseHeader(ifile)
    assert(module_list==['matplotlib', 'shilofue.Plot'])
    assert(object_list==[['pyplot'], ['LINEARPLOT']])


def test_scripting():
    '''
    Test the class SCRIPTING
    Assert:
        file is generated
    '''
    # test 1: old file
    ifile = os.path.join(source_dir, "PlotStatistics.py")
    ofile = os.path.join(test_dir, "PlotStatistics_scripting.py")
    assert(os.path.isfile(ifile))
    Scripting = SCRIPTING(ifile)
    Scripting(ofile)
    assert(os.path.isfile(ofile))
    # test 2: new file, check for file contents
    ifile = os.path.join(source_dir, "PlotStatistics1.py")
    ofile = os.path.join(test_dir, "PlotStatistics_scripting1.py")
    ofile_std = os.path.join(source_dir, "PlotStatistics1_std.py")
    assert(os.path.isfile(ifile))
    Scripting = SCRIPTING(ifile)
    Scripting(ofile)
    assert(os.path.isfile(ofile))
    assert(filecmp.cmp(ofile, ofile_std))  # assert file contents


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

