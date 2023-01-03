# -*- coding: utf-8 -*-
r"""Test for VtkPp.py, using big data files

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
import shilofue.VtkPp as VtkPp
from shilofue.TwoDSubduction0.VtkPp import *  # import test module
from shilofue.ParsePrm import ParseFromDealiiInput
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories
import vtk

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'TwoDSubduction', 'test_vtk_pp_slab')
TwoDSubduction_DIR = os.environ['TwoDSubduction_DIR']
has_project_root = (os.path.isdir(TwoDSubduction_DIR))
ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA = True


def test_slab_analysis():
    ''' 
    test the SlabAnalysis class, only works if the project files are presented
    '''
    source_dir1 = os.path.join(source_dir, 'slab_analysis') 
    # test 1: compute dynamic pressure by post-processing
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_std')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 105 # 10 Ma
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        SlabAnalysis(case_dir, vtu_snapshot, o_file, output_slab=True)
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))
    # test 2: with dynamic pressure outputed in vtu file.
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_dp_std')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces_dp")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 105 # 10 Ma
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3_dp', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        SlabAnalysis(case_dir, vtu_snapshot, o_file, output_slab=True)
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces_dp.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))
    # test 3: with dynamic pressure outputed in vtu file, use temperature differences
    # as the criteria for slab, dT = 100 K.
    # Note there are files containing more information generated, in "eba_cdpt_SA80.0_OA40.0/vtk_outputs"
    # e.g. the slab_env1_00101.vtp and slab_env0_00101.vtp contains a profile of the slab surface
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_dp_T100_std')
        o_env_std = os.path.join(source_dir1, 'slab_env0_00101_std.vtp')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces_dp_T100")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 106 # 10.1 Ma
        _time = 10.1e6
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3_dp', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        o_env = os.path.join(case_dir, "vtk_outputs", "slab_env0_00101.vtp")
        if os.path.isfile(o_env):
            os.remove(o_env)
        # get the interpolation function for a horizontally averaged T profile
        ha_file = os.path.join(case_dir, "output", "depth_average.txt")
        DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
        DepthAverage.Import(ha_file)
        itp_func = DepthAverage.GetInterpolateFunc(_time, "temperature")
        SlabAnalysis(case_dir, vtu_snapshot, o_file, use_dT=True, output_slab=True, output_poly_data=True, slab_envelop_interval=20e3)  # use the vertical profile from field data
        assert(os.path.isfile(o_env))
        assert(filecmp.cmp(o_env, o_env_std))  # compare file contents
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces_dp_T100.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))