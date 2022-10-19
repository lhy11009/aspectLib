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
from shilofue.TwoDSubduction0.VtkPp import get_theta
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
        self.slab_envelop_interval_z = kwargs.get("slab_envelop_interval_z", 10e3)
        self.slab_envelop_interval_y = kwargs.get("slab_envelop_interval_y", 10e3)
        self.slab_envelop_point_list0 = []
        self.slab_envelop_point_list1 = []
        self.slab_y_max = 0.0
        self.trench_coords_x = []
        self.trench_coords_y = []
        self.trench_coords_z = []

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
                if y > self.slab_y_max:
                    self.slab_y_max = y
        self.slab_depth = self.Ro - min_r  # cart
        print("PrepareSlabByPoints: %d points found in the subducting slab" % len(self.slab_points))
        # prepare the slab internal
        total_en_interval_z = int((self.slab_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval_z + 1)
        total_en_interval_y = int((self.slab_y_max) // self.slab_envelop_interval_y + 1)
        slab_en_point_lists = [ [] for i in range(total_en_interval_z * total_en_interval_y) ]
        for id in self.slab_points:
            x = points[id][0] # first, separate cells into intervals
            y = points[id][1]
            z = points[id][2]
            r = VtkPp.get_r3(x, y, z, self.geometry)
            id_en_z =  int(np.floor(
                                  (self.Ro - r - self.slab_shallow_cutoff)/
                                  self.slab_envelop_interval_z))# id in the envelop list
            id_en_y =  int(np.floor(y / self.slab_envelop_interval_y))
            id_en = id_en_y * total_en_interval_z + id_en_z
            slab_en_point_lists[id_en].append(id)
        trench_coords_x = []
        trench_coords_y = []
        trench_coords_z = []
        for id_en in range(len(slab_en_point_lists)):
            theta_min = 0.0  # then, loop again by intervals to look for a
            theta_max = 0.0  # max theta and a min theta for each interval
            point_list = slab_en_point_lists[id_en]
            if len(point_list) == 0:
                continue  # make sure we have some point
            is_first = True
            id_min = -1
            id_max = -1
            for id in point_list:
                x = points[id][0]
                y = points[id][1]
                theta = get_theta(x, y, self.geometry)  # cart
                if is_first:
                    id_min = id
                    id_max = id
                    theta_min = theta
                    theta_max = theta
                    is_first = False
                else:
                    if theta < theta_min:
                        id_min = id
                        theta_min = theta
                    if theta > theta_max:
                        id_max = id
                        theta_max = theta
            self.slab_envelop_point_list0.append(id_min)  # first half of the envelop
            self.slab_envelop_point_list1.append(id_max)  # second half of the envelop
            if id_en % total_en_interval_z == 0:
                # record trench coordinates
                id_en_y = id_en // total_en_interval_z
                if id_max > 0:
                    trench_coords_x.append(points[id_max][0])
                    trench_coords_y.append(points[id_max][1])
                    trench_coords_z.append(points[id_max][2])
        self.trench_coords_x = trench_coords_x
        self.trench_coords_y = trench_coords_y
        self.trench_coords_z = trench_coords_z
    
    def ExportSlabInfo(self):
        '''
        Output slab information
        '''
        return self.trench_coords_x, self.trench_coords_y, self.trench_coords_z


def SlabMorphology(case_dir, vtu_snapshot, **kwargs):
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
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    print("geometry: %s, Ro: %f, Xmax: %f" % (geometry, Ro, Xmax))
    # Initiate the working class
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax, slab_shallow_cutoff=70e3)
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
    # export slab points
    # 1. slab intervals
    fileout = os.path.join(output_path, 'slab.vtu')
    slab_point_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data, VtkP.slab_points)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_point_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab internal points: %s" % fileout)
    # 2. slab envelops
    fileout = os.path.join(output_path, 'slab_env0.vtu')
    slab_env0_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data,VtkP.slab_envelop_point_list0)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_env0_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab envelop0 (smaller x) points: %s" % fileout)
    fileout = os.path.join(output_path, 'slab_env1.vtu')
    slab_env1_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data,VtkP.slab_envelop_point_list1)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_env1_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab envelop1 (bigger x) points: %s" % fileout)
    end = time.time()
    print("Write slab points, takes %.2f s" % (end - start))
    start = end
    # prepare outputs
    trench_coords_x, trench_coords_y, trench_coords_z = VtkP.ExportSlabInfo()
    outputs = "# trench coordinates: x, y and z\n" + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(trench_coords_x)):
        outputs += str(trench_coords_x[i]) + " " + str(trench_coords_y[i]) + " " + str(trench_coords_z[i]) + "\n"
    fileout = os.path.join(output_path, "trench_%05d.txt" % vtu_snapshot)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions: %s" % fileout)
    end = time.time()
    print("Write trench positions, takes %.2f s" % (end - start))
    start = end

    # outputs = "%-12s%-12d%-14.4e\n"\
    # % (vtu_step, step, _time)
    # print("Output file generated: ", outputs)
    return vtu_step, outputs


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
        SlabMorphology(arg.inputs, int(arg.vtu_snapshot), rewrite=1)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()