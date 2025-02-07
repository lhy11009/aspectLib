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
# import filecmp  # for compare file contents
import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.PlotDepthAverage import *
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_depth_average')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_buoyancy():
    '''
    test function ComputeBuoyancy
    Asserts:
        plot utilities work
        values of the buoyancy calculated
    '''
    da_file0 = os.path.join(source_dir, "depth_average_1573.txt")
    da_file1 = os.path.join(source_dir, "depth_average_673.txt")
    fig_path = os.path.join(test_dir, "buoyancy.png")
    fig = plt.figure(tight_layout=True, figsize=(6, 10))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    if os.path.isfile(fig_path):
        os.remove(fig_path)

    # plot figures
    _, buoyancies, buoyancy_ratios = ComputeBuoyancy(da_file0, da_file1, ax1=ax1, ax2=ax2)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    # save figures
    fig_path = os.path.join(fig_path)
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))

    # check on the buoyancy calculated
    buoyancy440 = -562.7024359861579
    buoyancy_ratio440 = -0.16577675940292663
    assert(abs(buoyancies[440] - buoyancy440)/buoyancy440 < 1e-6)
    assert(abs(buoyancy_ratios[440] - buoyancy_ratio440)/buoyancy_ratio440 < 1e-6)
    pass

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

