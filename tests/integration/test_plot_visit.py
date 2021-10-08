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


# import os
# import pytest
import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
# from shilofue.Utilities import 
from shilofue.PlotVisit import *
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_visit')
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_visit')
test_cases_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_get_snaps_steps():
    case_dir = os.path.join(test_cases_dir, 'foo')
    
    # call function for graphical outputs
    snaps, times, steps = GetSnapsSteps(case_dir)
    # assertions
    assert(snaps == [6, 7, 8, 9])
    assert(times == [0.0, 100000.0, 200000.0, 300000.0])
    assert(steps == [0, 104, 231, 373])
    
    # call function for particle outputs
    snaps, times, steps = GetSnapsSteps(case_dir, 'particle')
    # assertions
    assert(snaps == [0, 1])
    assert(times == [0.0, 2e5])
    assert(steps == [0, 231])


def test_visit_options(): 
    # check visit_options (interpret script from standard ones)
    case_dir = os.path.join(test_cases_dir, 'test_visit')
    Visit_Options = VISIT_OPTIONS(case_dir)
    # call function
    Visit_Options.Interpret()
    ofile = os.path.join(test_dir, 'temperature.py')
    visit_script = os.path.join(source_dir, 'temperature.py')
    visit_script_base = os.path.join(source_dir, 'base.py')
    Visit_Options.read_contents(visit_script_base, visit_script)
    # make a new directory
    img_dir = os.path.join(test_dir, 'img')
    if os.path.isdir(img_dir):
        rmtree(img_dir)
    os.mkdir(img_dir)
    Visit_Options.options["IMG_OUTPUT_DIR"] = img_dir
    Visit_Options.substitute()
    ofile_path = Visit_Options.save(ofile)
    # assert file generated
    assert(os.path.isfile(ofile_path))
    # assert file is identical with standard
    ofile_std = os.path.join(source_dir, 'temperature_std.py')
    assert(os.path.isfile(ofile_std))
    assert(filecmp.cmp(ofile_path, ofile_std))

def test_vtk_TwoDSubduction_SlabAnalysis_options():
    # check options for vtk
    case_dir = os.path.join(test_cases_dir, 'test_vtk')
    option_path = os.path.join(test_dir, 'TwoDSubduction_SlabAnalysis.input')
    vtk_option_path = PrepareVTKOptions(case_dir, 'TwoDSubduction_SlabAnalysis', step=0, output=option_path)
    RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)


#def test_visit_options():
#    case_dir = os.path.join(test_source_dir, 'foo')
#    case_dir = os.path.join(test_cases_dir, 'foo')
#
#    # initiate 
#    Visit_Options = VISIT_OPTIONS(case_dir)
#
#    # call the Interpret function
#    Visit_Options.Interpret()
#
#    # compare to standard
#    json_file = os.path.join(case_dir, 'odict_std.json')
#    with open(json_file, 'r') as fin:
#        odict_std = json.load(fin)
#
#    assert(Visit_Options.odict == odict_std)
#    pass

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

