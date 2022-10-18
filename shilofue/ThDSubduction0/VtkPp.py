# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

descriptions

    This script handles analysis using the vtk package for this project
"""
#### 3rd parties 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from numpy import linalg as LA 
import multiprocessing
import time
#### self
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.ThDSubduction0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile
from shilofue.Plot import LINEARPLOT
import shilofue.VtkPp as VtkPp
# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
    analyze slab morphology, by parsing a pvtu file:\n\
        python -m shilofue.ThDSubduction0.VtkPp analyze_slab -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/gmg_test_stampede/test_ThD_gmg_mv1e20_picard_correction_side_plate -vs 0 \n\
\n\
")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
    '''
    def __init__(self, **kwargs):
        VtkPp.VTKP.__init__(self, **kwargs)  # initiation of parental class
        self.slab_shallow_cutoff = kwargs.get('slab_shallow_cutoff', 100e3)  # depth limit to slab
        self.slab_points = []  # points in the slab

    def PrepareSlabByPoints(self, slab_field_names, **kwargs):
        '''
        prepare slab composition, with position of points. I do it this way
        because the cell center data requires interpolation, and that takes about
        few minutes when my case scales up.
        '''
        slab_threshold = kwargs.get('slab_threshold', 0.2)
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        point_data = self.i_poly_data.GetPointData()
        # slab composition field
        slab_field = VtkPp.OperateDataArrays(point_data, slab_field_names,\
        [0 for i in range(len(slab_field_names) - 1)])
        # add cells by composition
        min_r = self.Ro
        for i in range(self.i_poly_data.GetNumberOfPoints()):
            x = points[i][0]
            y = points[i][1]
            z = points[i][2]
            r = VtkPp.get_r3(x, y, z, self.geometry)
            slab = slab_field[i]
            if slab > slab_threshold and ((self.Ro - r) > self.slab_shallow_cutoff):
                self.slab_points.append(i)
                if r < min_r:
                    min_r = r
        self.slab_depth = self.Ro - min_r  # cart
        # todo
        # export the internal points

    def ExportSlabInternalPointsToFile(self, fileout):
        '''
        Export the coordinates of slab internal points to file
        '''
        slab_point_grid = self.ExportSlabInternalPoints()
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputData(slab_point_grid)
        writer.SetFileName(fileout)
        writer.Update()
        writer.Write()
        print("ExportSlabInternalPoints: write file %s" % fileout)
    
    def ExportSlabInternalPoints(self, output_xy=False):
        '''
        export slab internal points
        '''
        assert(self.slab_points is not [])
        slab_vtk_points = vtk.vtkPoints()
        for id in self.slab_points:
            xs = self.i_poly_data.GetPoint(id)
            slab_vtk_points.InsertNextPoint(xs[0], xs[1], xs[2])
        slab_point_grid = vtk.vtkUnstructuredGrid()
        slab_point_grid.SetPoints(slab_vtk_points)
        if output_xy:
            coords = vtk_to_numpy(slab_point_grid.GetPoints().GetData())
            return coords
        else:
            return slab_point_grid


def SlabAnalysis(case_dir, vtu_snapshot, **kwargs):
    '''
    Wrapper for using PVTK class to analyze one step
    Inputs:
        case_dir (str): case directory
        vtu_step (int): step in vtu outputs
    '''
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    assert(os.path.isfile(filein))
    # output_path
    output_path = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    # Initiate an object for the options we use for vtk, also get time and step
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    # Initiate the working class
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    # call some functions
    field_names = ['T', 'p', 'density', 'sp_upper', 'sp_lower']
    start = time.time() # record time
    VtkP.ConstructPolyData(field_names, include_cell_center=False)
    end = time.time()
    print("Construct polydata, takes %.2f s" % (end - start))
    start = end
    # Analyze slab composiiotn by points
    VtkP.PrepareSlabByPoints(['sp_upper', 'sp_lower'])
    end = time.time()
    print("Prepare slab composition, takes %.2f s" % (end - start))
    start = end
    # output the slab internal points
    fileout = os.path.join(output_path, 'slab.vtu')
    VtkP.ExportSlabInternalPointsToFile(fileout)
    end = time.time()
    print("Write slab points, takes %.2f s" % (end - start))
    start = end
    # generate some outputs
    # outputs = "%-12s%-12d%-14.4e\n"\
    # % (vtu_step, step, _time)
    # print("Output file generated: ", outputs)
    # return vtu_step, outputs


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='Some outputs')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='step')
    parser.add_argument('-vs', '--vtu_step', type=int,
                        default=0,
                        help='vtu_step')
    parser.add_argument('-vss', '--vtu_snapshot', type=int,
                        default=0,
                        help='vtu_snapshot')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'analyze_slab':
        # use the wrapper
        SlabAnalysis(arg.inputs, int(arg.vtu_snapshot), rewrite=1)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()