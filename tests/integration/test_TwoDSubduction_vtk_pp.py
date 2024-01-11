# -*- coding: utf-8 -*-
r"""Test for VtkPp.py

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


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_GetTimeDepthTip():
    '''
    test GetTimeDepthTip
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp', "test_TwoD_vtk_case1")
    SlabPlot = SLABPLOT('test')
    t_660_std = 2419872.042534338  # yr
    t_660 = SlabPlot.GetTimeDepthTip(case_dir, 660e3)
    assert(abs(t_660 - t_660_std) / t_660_std < 1e-6)


def test_wedge_T():
    '''
    Test wedge temperature
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_wedge_T")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # test 1 output slab grid & envelop
    fileout = os.path.join(output_path, 'wedge_T100.txt')
    VtkP.ExportWedgeT(fileout=fileout)
    fileout_std = os.path.join(case_dir, 'wedge_T100_std.txt')
    assert(os.path.isfile(fileout))
    assert(filecmp.cmp(fileout_std, fileout))  # compare file extent



def test_prepare_slab():
    '''
    Test utilities from VtkPp.py
    0: read in pvtu file
    test 1: output slab grid
    '''
    # assert something 
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_readfile")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # test 1 output slab grid & envelop
    fileout = os.path.join(output_path, 'slab.vtu')
    fileout_std = os.path.join(case_dir, 'slab_std.vtu')
    slab_grid = VtkP.ExportSlabInternal()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    assert(os.path.isfile(fileout))  # assert file existence
    assert(filecmp.cmp(fileout_std, fileout))  # compare file extent
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord() # envelop
    assert(abs(slab_envelop0[0, 0] - 5086593.875)/5086593.875 < 1e-8)
    assert(abs(slab_envelop0[6, 1] - 3723199.125)/3723199.125< 1e-8)
    # test 2, slab buoyancy
    r0_range = [6371e3 - 2890e3, 6371e3]
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    total_buoyancy, b_profile = VtkP.SlabBuoyancy(v_profile, 5e3)
    assert(abs(total_buoyancy + 5394703810473.24)/5394703810473.24 < 1e-8)
    assert((b_profile[11, 1]-1.09996017e+12)/1.09996017e+12 < 1e-8)


def test_export_slab_info_sph():
    '''
    Test slab properties from VtkPp.py
    assert:
        1. value of trench position, slab_depth, dip angle
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_readfile")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    assert(abs(trench - 0.6466922210397674)/0.6466922210397674 < 1e-6)
    assert(abs(slab_depth - 191927.42159304488)/191927.42159304488 < 1e-6)
    assert(abs(dip_100 - 1.061791071552635)/1.061791071552635< 1e-6)


def test_export_slab_info_cart():
    '''
    Test slab properties from VtkPp.py
    assert:
        1. value of trench position, slab_depth, dip angle
    Note:
        here the error in dip100 has a lot to do with resolution. If
        we take both the 5th adaptive refinement, both of these angles (chunk and box)
        are around 0.63
    '''
    case_dir = os.path.join(source_dir, 'cartesian')
    prm_file = os.path.join(case_dir, 'case.prm')
    assert(os.path.isfile(prm_file))
    with open(prm_file, 'r') as fin:
        idict = ParseFromDealiiInput(fin)
    geometry = idict['Geometry model']['Model name']
    if geometry == 'chunk':
        Ro = float(idict['Geometry model']['Chunk']['Chunk outer radius'])
    elif geometry == 'box':
        Ro = float(idict['Geometry model']['Box']['Y extent'])
    output_path = os.path.join(test_dir, "TwoDSubduction_vtk_pp_slab")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP(geometry=geometry, Ro=Ro)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    assert(abs(trench - 4210541.75)/4210541.75 < 1e-6)
    assert(abs(slab_depth - 220136.75)/220136.75 < 1e-6)
    assert(abs(dip_100- 0.6702772823940486)/0.6702772823940486 < 1e-6)


def test_export_velocity():
    '''
    test function ExportVelocity
    assert:
        1. velocity of the overiding plate and the subducting plate
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_readfile")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    vsp, vov = VtkP.ExportVelocity()
    assert(abs(vsp[0]) < 1e-6 and abs(vsp[1]) < 1e-6 and abs(vsp[2]) < 1e-6)  # thest two values are 0.0
    assert(abs(vov[0]) < 1e-6 and abs(vov[1]) < 1e-6 and abs(vov[2]) < 1e-6)
    # assert


def test_analysize_sz():
    '''
    test the sz geometry
    '''
    case_dir = os.path.join(ASPECT_LAB_DIR, "tests/integration/fixtures/big_files/test_TwoD_vtk_pp_full")
    assert(os.path.isdir(case_dir))
    # fix the output directory
    output_path = os.path.join(test_dir, "twod_vtkp_sz")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    vtu_snapshot = 97
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("SlabMorphology: processing %s" % filein)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    # vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    # call the functions for the shear zone
    fileout = os.path.join(output_path, "sz.txt")
    file_std = os.path.join(case_dir, "sz_std.txt")
    VtkP.PrepareSZ(fileout)
    assert(os.path.isfile(fileout))  # assert file generation
    assert(filecmp.cmp(fileout, file_std))  # compare file contents
    # plot
    fig_path = os.path.join(output_path, "sz_thickness.png") 
    fig, ax = plt.subplots()
    MorphPlotter = SLABPLOT("plot_slab")
    MorphPlotter.PlotShearZoneThickness(case_dir, trench, axis=ax, filein=fileout, label='shear zone thickness')
    ax.legend()
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))  # assert figure generation


def test_find_mdd():
    case_dir = os.path.join(ASPECT_LAB_DIR, "tests/integration/fixtures/big_files/test_TwoD_vtk_pp_full")
    assert(os.path.isdir(case_dir))
    vtu_snapshot = 97
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("SlabMorphology: processing %s" % filein)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    # vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    mdd = VtkP.FindMDD()
    assert(abs(mdd - 7.9729e+04) / 7.9729e+04 < 1e-3) # check the mdd value
    pass


def test_trench_T():
    '''
    Test trench temperature
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_trench_T")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # test 1 output slab grid & envelop
    fileout = os.path.join(output_path, 'trench_T.txt')
    VtkP.ExportTrenchT(fileout=fileout)
    fileout_std = os.path.join(case_dir, 'trench_T_std.txt')
    assert(os.path.isfile(fileout))
    assert(filecmp.cmp(fileout_std, fileout))  # compare file extent


def test_fit_trench_T():
    '''
    test fitting an age of the subducting plate to the trench temperature profile
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(case_dir, 'vtk_outputs')
    if os.path.isdir(output_path):
        rmtree(output_path)
    vtu_snapshot = 2
    # export the trench temperature
    TrenchT(case_dir, vtu_snapshot)
    SlabPlot = SLABPLOT('trench_T')
    age_myr_std = 77.78356193082101
    age_myr = SlabPlot.FitTrenchT(case_dir, vtu_snapshot)
    assert(abs(age_myr - age_myr_std)/age_myr_std < 1e-6)
