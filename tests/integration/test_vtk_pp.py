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
from shilofue.VtkPp import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')



if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_utilities():
    '''
    Test utilities from VtkPp.py
    0: read in pvtu file
    1: export to a ascii file
    2: extract a temperature contour
    3: interpolation
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
    field_names = ['T', 'density']
    VtkP.ConstructPolyData(field_names)
    # test 1: output to ascii file
    fileout = os.path.join(output_path, 'ascii_output.txt')
    fileout_std = os.path.join(case_dir, 'ascii_output_std.txt')
    ExportPolyDataAscii(VtkP.i_poly_data, ['T', 'density'], fileout)
    assert(os.path.isfile(fileout))
    # assert the contents of file
    assert(filecmp.cmp(fileout, fileout_std))
    # test 2: extract temperature contour at 1200 C (1473 K)
    fileout = os.path.join(output_path, 'contour_output.txt')
    fileout_std = os.path.join(case_dir, 'contour_output_std.txt')
    ExportContour(VtkP.i_poly_data, 'T', 1000.0, fileout=fileout)
    assert(os.path.isfile(fileout))
    # assert the contents of file
    assert(filecmp.cmp(fileout, fileout_std))
    # test 3: fields at theta = 1 deg
    n = 10
    ro = 6371e3  # earth's outer radius
    ri = 6371e3 - 2890e3  # earth's inner radius
    thetas = np.ones(n) * np.pi / 180.0
    rs = np.linspace(ri + 1e3, ro - 1e3, n)
    xs = rs * np.cos(thetas)
    ys = rs * np.sin(thetas)
    ps = np.zeros((n, 2))
    ps[:, 0] = xs
    ps[:, 1] = ys
    o_poly_data = InterpolateGrid(VtkP.i_poly_data, ps)
    o_point_data = o_poly_data.GetPointData()
    densities = vtk_to_numpy(o_point_data.GetArray('density'))
    Ts = vtk_to_numpy(o_point_data.GetArray('T'))
    assert(abs(densities[0] - 3761.43)/3761.43 < 1e-6)  # assert field data
    assert(abs(densities[1] - 3772.0)/3772.0 < 1e-6)
    assert(abs(Ts[0] - 3493.1062) / 3493.1062 < 1e-6)
    assert(abs(Ts[1] - 2971.1995) / 2971.1995 < 1e-6)


def test_static_pressure():
    '''
    test the function of StaticPressure
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density']
    VtkP.ConstructPolyData(field_names)
    # extract static pressure
    ro = 6371e3
    r = ro - 100e3
    theta1 = 1  # under the subducting plate
    theta2 = 61 * np.pi / 180.0  # under the overiding plate
    static_pressure_std1 = 3.228537216558201e9
    static_pressure1 = VtkP.StaticPressure([r, ro-10.0], theta1, 500, grav_acc=9.8)  # n = 500 - 1000
    assert((static_pressure1-static_pressure_std1)/ static_pressure_std1 <  1e-6)
    static_pressure_std2 = 3.224576241711115e9 # smaller than previous, as the subducting plate is older.
    static_pressure2 = VtkP.StaticPressure([r, ro-10.0], theta2, 500, grav_acc=9.8)  # n = 500 - 1000
    assert((static_pressure2-static_pressure_std2)/ static_pressure_std2 < 1e-6)
    # assert(static_pressure == )


def test_gravity_data():
    '''
    test the function of importing gravity
    Assert:
        1. value of gravity at 1000 km depth
    '''
    aspect_source_dir = Utilities.var_subs('${ASPECT_SOURCE_DIR}')
    gravity_file = os.path.join(aspect_source_dir, "data", "gravity-model", "prem.txt")  # default gravity file (prem) for aspect
    assert(os.path.isfile(gravity_file))
    VtkP = VTKP()
    VtkP.ImportGravityData(gravity_file)
    grav_acc_1000 = VtkP.GetGravityAcc(1e6)
    assert((grav_acc_1000 - 9.966665)/9.966665 < 1e-6)


# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

