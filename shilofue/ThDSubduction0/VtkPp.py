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
import math
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from matplotlib import pyplot as plt
from matplotlib import cm, gridspec
from numpy import linalg as LA 
import multiprocessing
import time
from joblib import Parallel, delayed
import warnings
from scipy.interpolate import interp1d, UnivariateSpline, griddata
#### self
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.ThDSubduction0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile
from shilofue.Plot import LINEARPLOT
from shilofue.PlotCombine import PLOT_COMBINE
from shilofue.TwoDSubduction0.VtkPp import get_theta
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
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
        trench morphology\n\
        python -m shilofue.ThDSubduction0.VtkPp morph_case -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/EBA_2d_consistent_8_6/eba3d_width51_c22_AR4 -ti 1e6 \n\
\n\
        run morphology in parallel\n\
        python -m shilofue.ThDSubduction0.VtkPp morph_case_parallel -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/EBA_2d_consistent_8_6/eba3d_width51_c22_AR4 -ti 1e6 \n\
\n\
    plot trench positions, assuming the trench.txt files have already been generated: \n\
        the -ti option gives a time interval\n\
        python -m shilofue.ThDSubduction0.VtkPp plot_trench -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/gmg_test_stampede/test_ThD_gmg_mv1e20_picard_correction_side_plate -ti 2e6 \n\
\n\
        plot trench morphology in episodes\n\
        python -m shilofue.ThDSubduction0.VtkPp  plot_trench_by_episodes -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/EBA_2d_consistent_8_6/eba3d_width51_c22_AR4\n\
\n\
")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
    '''
    def __init__(self, **kwargs):
        VtkPp.VTKP.__init__(self, **kwargs)  # initiation of parental class
        self.slab_shallow_cutoff = kwargs.get('slab_shallow_cutoff', 70e3)  # depth limit to slab
        self.slab_points = []  # points in the slab
        self.slab_envelop_interval_z = kwargs.get("slab_envelop_interval_z", 10e3)
        self.slab_envelop_interval_y = kwargs.get("slab_envelop_interval_y", 10e3)
        self.slab_envelop_point_list0 = []
        self.slab_envelop_point_list1 = []
        self.slab_y_max = 0.0
        self.trench_coords_x = []
        self.trench_coords_y = []
        self.trench_coords_z = []
        self.trench_b_coords_x = []
        self.trench_b_coords_y = []
        self.trench_b_coords_z = []
        self.center_profile_x = []
        self.center_profile_y = []
        self.center_profile_z = []
        self.pinD_profile_x = []
        self.pinD_profile_y = []
        self.pinD_profile_z = []
        self.pinD_b_profile_x = []
        self.pinD_b_profile_y = []
        self.pinD_b_profile_z = []
        self.slab_depth = 0.0  # slab depth
        self.buoyancies = None

    def PrepareSlabByPoints(self, slab_field_names, **kwargs):
        '''
        prepare slab composition, with position of points. I do it this way
        because the cell center data requires interpolation, and that takes about
        few minutes when my case scales up.
        '''
        slab_threshold = kwargs.get('slab_threshold', 0.2)
        field_type = kwargs.get('field_type', "composition")
        # todo_pind
        pin_depth = kwargs.get("pin_depth", 100e3)
        get_pin_depth = kwargs.get("get_pin_depth", 0)
        if field_type == "composition":
            # tranferred to to enum type to save computational cost
            field_enum = 0
        elif field_type == "temperature":
            field_enum = 1
        else:
            raise ValueError("field_type needs to be either composition or temperature")
        if get_pin_depth:
            # the pin depth must be deeper than the cutoff depth to ensure successful pinning
            assert(pin_depth > self.slab_shallow_cutoff)
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
            if field_enum == 0:
                # with this method, I use the total value of a few composition fields
                # and determine slab internal points based on higher compositions.
                condition = (slab > slab_threshold)
            elif field_enum == 1:
                # with this method, I use the value of the temperature field
                # and determine slab internal points based on lower temperature.
                condition = (slab < slab_threshold)
            if condition and ((self.Ro - r) > self.slab_shallow_cutoff):
                # there is also a shallow cutoff applied.
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
        id_en_pin_depth = int((pin_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval_z)
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
        trench_b_coords_x = []
        trench_b_coords_y = []
        trench_b_coords_z = []
        center_profile_x = []
        center_profile_y = []
        center_profile_z = []
        pinD_coords_x = []
        pinD_coords_y = []
        pinD_coords_z = []
        pinD_b_coords_x = []
        pinD_b_coords_y = []
        pinD_b_coords_z = []
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
            if id_en < total_en_interval_z:
                # record center profile
                # todo_center
                if id_max > 0:
                    center_profile_x.append(points[id_max][0])
                    center_profile_y.append(points[id_max][1])
                    center_profile_z.append(points[id_max][2])
            if id_en % total_en_interval_z == 0:
                # record trench coordinates
                id_en_y = id_en // total_en_interval_z
                if id_max > 0:
                    trench_coords_x.append(points[id_max][0])
                    trench_coords_y.append(points[id_max][1])
                    trench_coords_z.append(points[id_max][2])
                if id_min > 0:
                    trench_b_coords_x.append(points[id_min][0])
                    trench_b_coords_y.append(points[id_min][1])
                    trench_b_coords_z.append(points[id_min][2])
            if get_pin_depth:
                if id_en % total_en_interval_z == id_en_pin_depth:
                    # get the pinned x, y and z at a pinned depth
                    if id_max > 0:
                        pinD_coords_x.append(points[id_max][0])
                        pinD_coords_y.append(points[id_max][1])
                        pinD_coords_z.append(points[id_max][2])
                    if id_min > 0:
                        pinD_b_coords_x.append(points[id_min][0])
                        pinD_b_coords_y.append(points[id_min][1])
                        pinD_b_coords_z.append(points[id_min][2])
        # append data to class objects
        self.trench_coords_x = trench_coords_x
        self.trench_coords_y = trench_coords_y
        self.trench_coords_z = trench_coords_z
        self.trench_b_coords_x = trench_b_coords_x
        self.trench_b_coords_y = trench_b_coords_y
        self.trench_b_coords_z = trench_b_coords_z
        self.center_profile_x = center_profile_x
        self.center_profile_y = center_profile_y
        self.center_profile_z = center_profile_z
        self.pinD_profile_x = pinD_coords_x
        self.pinD_profile_y = pinD_coords_y
        self.pinD_profile_z = pinD_coords_z
        self.pinD_b_profile_x = pinD_b_coords_x
        self.pinD_b_profile_y = pinD_b_coords_y
        self.pinD_b_profile_z = pinD_b_coords_z

    # todo_hv 
    def SlabSurfaceAtDepth(self, depth):
        '''
        Get the slab surface at a given depth
        Here I get the slab surface in a cross section at a given depth
        Inputs:
            depth - a given depth to get a cross section
        '''
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        assert(len(self.slab_envelop_point_list0) > 0)  # check the slab envelops have points and depth within the range
        assert(depth < self.slab_depth)
        xs = []
        ys = []
        zs = []
        # here I don't proceed with the exact id depth, 
        # because there are missing intervals in the slab_envelop_point_list0
        # (otherwise determine a i_z -> figure out the i that way)
        for i in range(len(self.slab_envelop_point_list0)):
            id = self.slab_envelop_point_list0[i]
            x = points[id][0]
            y = points[id][1]
            z = points[id][2]
            if abs(self.Ro - depth - z) < 10e3:
                xs.append(x)
                ys.append(y)
                zs.append(z)
        return xs, ys, zs


    def ComputeBuoyancy(self):
        '''
        Compute buoyancy
        '''
        # initiate
        grav_acc = 10.0
        density_data = vtk_to_numpy(self.i_poly_data.GetPointData().GetArray('density'))
        self.buoyancies = np.zeros(len(self.slab_points))
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())

        # compute buoyancy
        for i in range(len(self.slab_points)):
            id_en = self.slab_points[i]
            x = points[id_en][0] # first, separate cells into intervals
            y = points[id_en][1]
            z = points[id_en][2]
            r = VtkPp.get_r3(x, y, z, self.geometry)
            depth = self.Ro - r
            density = density_data[id_en]
            density_ref = self.density_ref_func(depth)
            self.buoyancies[i] = - grav_acc * (density - density_ref)  # gravity
        pass
    
    def ExportSlabInfo(self):
        '''
        Output slab information
        '''
        return self.trench_coords_x, self.trench_coords_y, self.trench_coords_z

    def ExportSlabInfoB(self):
        '''
        Output slab information (bottom)
        '''
        return self.trench_b_coords_x, self.trench_b_coords_y, self.trench_b_coords_z
    
    def ExportPinDInfo(self):
        '''
        Output pinned depth slab curve
        '''
        return self.pinD_profile_x, self.pinD_profile_y, self.pinD_profile_z
    
    def ExportPinDInfoB(self):
        '''
        Output pinned depth slab curve (bottom)
        '''
        return self.pinD_b_profile_x, self.pinD_b_profile_y, self.pinD_b_profile_z

    # todo_center 
    def ExportCenterProfile(self):
        '''
        Output slab information
        '''
        return self.center_profile_x, self.center_profile_y, self.center_profile_z


def SlabMorphology(case_dir, vtu_snapshot, **kwargs):
    '''
    Wrapper for using PVTK class to analyze one step
    Inputs:
        case_dir (str): case directory
        vtu_step (int): step in vtu outputs
    '''
    # inputs
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found" % filein)
    slab_envelop_interval_y = kwargs.get("slab_envelop_interval_y", 10e3)
    slab_envelop_interval_z = kwargs.get("slab_envelop_interval_z", 10e3)
    slab_shallow_cutoff = kwargs.get("slab_shallow_cutoff", 70e3)
    include_force_balance = kwargs.get("include_force_balance", False)
    crust_only = kwargs.get("crust_only", 0)
    horizontal_velocity_depths = kwargs.get("horizontal_velocity_depths", []) # provide a list of depth to investigate the velocity field
    pin_depth = kwargs.get("pin_depth", 100e3)
    appendix = "_%05d" % vtu_snapshot
    # output_path
    output_path = kwargs.get("output", os.path.join(case_dir, 'vtk_outputs'))
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
    if include_force_balance:
        ha_file = os.path.join(case_dir, "output", "depth_average.txt")
        assert(os.path.isfile(ha_file))
    else:
        ha_file = None

    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax, slab_shallow_cutoff=slab_shallow_cutoff,\
        slab_envelop_interval_y=slab_envelop_interval_y, slab_envelop_interval_z=slab_envelop_interval_z,\
        ha_file=ha_file, time=_time)
    VtkP.ReadFile(filein)
    
    # call some functions
    field_names = ['T', 'p', 'density', 'sp_upper', 'sp_lower']
    start = time.time() # record time
    VtkP.ConstructPolyData(field_names, include_cell_center=False, )
    end = time.time()
    print("Construct polydata, takes %.2f s" % (end - start))
    start = end
    # Analyze slab composiiotn by points
    # by default, we turn on the functionality to export slab curve at
    # a certain depth
    if crust_only == 1:
        slab_query_fields = ['sp_upper']
    elif crust_only == 2:
        slab_query_fields = ['sp_lower']
    else:
        slab_query_fields = ['sp_upper', 'sp_lower']
    VtkP.PrepareSlabByPoints(slab_query_fields, get_pin_depth=1, pin_depth=pin_depth)
    end = time.time()
    print("Prepare slab composition, takes %.2f s" % (end - start))
    start = end
    # output the slab internal points
    # export slab points
    # 1. slab intervals
    # fileout = os.path.join(output_path, 'slab.vtu')
    slab_point_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data, VtkP.slab_points)
    if include_force_balance:
        # mark the composition used for the envelop in the filename
        if crust_only == 1:
            filename = 'slab.vtp'
        elif crust_only == 2:
            filename = 'slab_l.vtp'
        else:
            filename = 'slab_lu.vtp'
        fileout = os.path.join(output_path, filename)
        # buoyancy array
        VtkP.ComputeBuoyancy()
        buoyancy_vtk_array = numpy_to_vtk(VtkP.buoyancies)
        buoyancy_vtk_array.SetName("buoyancy")
        o_poly_data = vtk.vtkPolyData()
        o_poly_data.SetPoints(slab_point_grid.GetPoints())
        o_poly_data.GetPointData().SetScalars(buoyancy_vtk_array)
        # density array
        densities = vtk_to_numpy(VtkP.i_poly_data.GetPointData().GetArray('density'))
        densities_slab = np.zeros(len(VtkP.slab_points))
        for i in range(len(VtkP.slab_points)):
            id_en = VtkP.slab_points[i]
            densities_slab[i] = densities[id_en]
        density_vtk_array = numpy_to_vtk(densities_slab)
        density_vtk_array.SetName("density")
        o_poly_data.GetPointData().AddArray(density_vtk_array)
        # horiz_average array, debug
        points = vtk_to_numpy(VtkP.i_poly_data.GetPoints().GetData())
        densities_ref = np.zeros(len(VtkP.slab_points))
        for i in range(len(VtkP.slab_points)):
            id_en = VtkP.slab_points[i]
            x = points[id_en][0] # first, separate cells into intervals
            y = points[id_en][1]
            z = points[id_en][2]
            r = VtkPp.get_r3(x, y, z, VtkP.geometry)
            depth = Ro - r
            densities_ref[i] = VtkP.density_ref_func(depth)
        densities_ref_vtk_array = numpy_to_vtk(densities_ref)
        densities_ref_vtk_array.SetName("density_ref")
        o_poly_data.GetPointData().AddArray(densities_ref_vtk_array)
        # export data
        VtkPp.ExportPolyData(o_poly_data, fileout)
    else:
        if crust_only == 1:
            filename = 'slab.vtu'
        elif crust_only == 2:
            filename = 'slab_l.vtu'
        else:
            filename = 'slab_lu.vtu'
        fileout = os.path.join(output_path, filename)
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputData(slab_point_grid)
        writer.SetFileName(fileout)
        writer.Update()
        writer.Write()
    print("File output for slab internal points: %s" % fileout)
    # 2. slab envelops
    # 2a bottom
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "slab_env0" + appendix + ".vtu"
    elif crust_only == 2:
        filename = "slab_env0_l" + appendix + ".vtu"
    else:
        filename = "slab_env0_lu" + appendix + ".vtu"
    fileout = os.path.join(output_path, filename)
    slab_env0_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data,VtkP.slab_envelop_point_list0)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_env0_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab envelop0 (smaller x) points: %s" % fileout)
    # 2b top
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "slab_env1" + appendix + ".vtu"
    elif crust_only == 2:
        filename = "slab_env1_l" + appendix + ".vtu"
    else:
        filename = "slab_env1_lu" + appendix + ".vtu"
    fileout = os.path.join(output_path, filename)
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
    # First write the position of the trench
    # a. top
    trench_coords_x, trench_coords_y, trench_coords_z = VtkP.ExportSlabInfo()
    outputs = "# trench coordinates: x, y and z\n" + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(trench_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(trench_coords_x[i]) + " " + str(trench_coords_y[i]) + " " + str(trench_coords_z[i])
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "trench_%05d.txt" % vtu_snapshot
    elif crust_only == 2:
        filename = "trench_l_%05d.txt" % vtu_snapshot
    else:
        filename = "trench_lu_%05d.txt" % vtu_snapshot
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions (upper boundary): %s" % fileout)
    # b. bottom
    trench_b_coords_x, trench_b_coords_y, trench_b_coords_z = VtkP.ExportSlabInfoB()
    outputs = "# trench coordinates: x, y and z\n" + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(trench_b_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(trench_b_coords_x[i]) + " " + str(trench_b_coords_y[i]) + " " + str(trench_b_coords_z[i])
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "trench_b_%05d.txt" % vtu_snapshot
    elif crust_only == 2:
        filename = "trench_l_b_%05d.txt" % vtu_snapshot
    else:
        filename = "trench_lu_b_%05d.txt" % vtu_snapshot
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)  
    print("File output for trench positions (lower boundary): %s" % fileout)
    end = time.time()
    print("Write trench positions, takes %.2f s" % (end - start))
    start = end
    # Then write the profile at the center of the slab
    center_profile_x, center_profile_y, center_profile_z = VtkP.ExportCenterProfile()
    outputs = "# center profile coordinates: x, y and z\n" + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(center_profile_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(center_profile_x[i]) + " " + str(center_profile_y[i]) + " " + str(center_profile_z[i])
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "center_profile_%05d.txt" % vtu_snapshot
    elif crust_only == 2:
        filename = "center_profile_l_%05d.txt" % vtu_snapshot
    else:
        filename = "center_profile_lu_%05d.txt" % vtu_snapshot
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for the center profile: %s" % fileout)
    end = time.time()
    print("Write center profile, takes %.2f s" % (end - start))
    start = end
    # Then write the pinned depth if required
    # top
    pinD_coords_x, pinD_coords_y, pinD_coords_z = VtkP.ExportPinDInfo()
    outputs = "# pinned points at depth %.2f km coordinates: x, y and z\n" % (pin_depth/1e3) + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(pinD_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(pinD_coords_x[i]) + " " + str(pinD_coords_y[i]) + " " + str(pinD_coords_z[i])
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "trench_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    elif crust_only == 2:
        filename = "trench_l_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    else:
        filename = "trench_lu_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions at depth %.2f km: %s" % (pin_depth/1e3, fileout))
    # bottom
    pinD_b_coords_x, pinD_b_coords_y, pinD_b_coords_z = VtkP.ExportPinDInfoB()
    outputs = "# pinned points at depth %.2f km coordinates: x, y and z\n" % (pin_depth/1e3) + \
    "# vtu step: %d\n" % vtu_step  + \
    "# time: %.4e\n" % _time
    for i in range(len(pinD_b_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(pinD_b_coords_x[i]) + " " + str(pinD_b_coords_y[i]) + " " + str(pinD_b_coords_z[i])
    if crust_only == 1:
        # mark the composition used for the envelop in the filename
        filename = "trench_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    elif crust_only == 2:
        filename = "trench_l_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    else:
        filename = "trench_lu_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions at depth %.2f km: %s" % (pin_depth/1e3, fileout))
    end = time.time()
    print("Write trench positions, takes %.2f s" % (end - start))
    start = end
    # Write a slice at a fixed depth
    # todo_hv
    for hv_depth in horizontal_velocity_depths:
        coords_x, coords_y, coords_z = VtkP.SlabSurfaceAtDepth(hv_depth)
        outputs = "# trench coordinates: x, y and z\n" + \
        "# vtu step: %d\n" % vtu_step  + \
        "# time: %.4e\n" % _time
        for i in range(len(coords_x)):
            if i > 0:
                outputs += "\n"
            outputs += str(coords_x[i]) + " " + str(coords_y[i]) + " " + str(coords_z[i])
        if crust_only == 1:
            # mark the composition used for the envelop in the filename
            filename = "slab_surface_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth/1e3)
        elif crust_only == 2:
            filename = "slab_surface_l_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth/1e3)
        else:
            filename = "slab_surface_lu_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth/1e3)
        fileout = os.path.join(output_path, filename)
        with open(fileout, 'w') as fout:
            fout.write(outputs)
        print("File output for trench positions: %s" % fileout)
    end = time.time()
    print("Write slab surface at depth, takes %.2f s" % (end - start))

    return vtu_step, ""
 
####
# Case-wise functions
####
def SlabMorphologyCase(case_dir, **kwargs):
    '''
    run vtk and get outputs for every snapshots
    Inputs:
        kwargs:
            rewrite: if rewrite previous results
    '''
    if_rewrite = kwargs.get('rewrite', 0)
    step_for_derivatives = kwargs.get('step_for_derivatives', 5)
    slab_envelop_interval_y = kwargs.get("slab_envelop_interval_y", 10e3)
    slab_envelop_interval_z = kwargs.get("slab_envelop_interval_z", 10e3)
    slab_shallow_cutoff = kwargs.get("slab_shallow_cutoff", 70e3)
    include_force_balance = kwargs.get("include_force_balance", False)
    use_parallel = kwargs.get('use_parallel', False)
    crust_only = kwargs.get("crust_only", False)
    # get all available snapshots
    # the interval is choosen so there is no high frequency noises
    time_interval_for_slab_morphology = kwargs.get("time_interval", 0.5e6)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
    # available_pvtu_steps = [i - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) for i in available_pvtu_snapshots]
    # get where previous session ends
    vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(vtk_output_dir):
        os.mkdir(vtk_output_dir)
    
    # slab_morph_file = os.path.join(vtk_output_dir, 'slab_morph.txt')

    # Initiation Wrapper class for parallel computation
    if crust_only == 1:
        # mark the composition used in the file name
        base_name = "slab_morph"
    elif crust_only == 2:
        base_name = "slab_morph_l"
    elif crust_only == 3:
        base_name = "slab_morph_lu"
    # use a wrapper of a processing function
    # the base_name is the basename of files to generate
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK(base_name, SlabMorphology, if_rewrite=if_rewrite, slab_envelop_interval_y=slab_envelop_interval_y,\
                                                slab_envelop_interval_z=slab_envelop_interval_z, slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only,\
                                                include_force_balance=include_force_balance)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file

#    if if_rewrite:
#        if os.path.isfile(slab_morph_file):
#            print("%s: Delete old slab_morph.txt file." % Utilities.func_name())
#            os.remove(slab_morph_file)  # delete slab morph file
#        ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite

    # loop for all the steps to plot
    # Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_step)\
    # for pvtu_step in available_pvtu_steps)  # first run in parallel and get stepwise output
    ParallelWrapper.clear()
    if use_parallel:
        num_cores = 2 # multiprocessing.cpu_count()
        print("run in parallel, number of cores: %d" % num_cores)
        # raise NotImplementedError("Parallel for the function %s is not properly implemented yet" % Utilities.func_name())
        Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_snapshot)\
        for pvtu_snapshot in available_pvtu_snapshots)  # first run in parallel and get stepwise output
        print("call assemble_parallel")  # debug
        pvtu_steps_o, outputs = ParallelWrapper.assemble_parallel()
    else:
        for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
            ParallelWrapper(pvtu_snapshot)
            # ParallelWrapper(pvtu_snapshot + step_for_derivatives)  # for computing the derivatives, might be bugs
        pvtu_steps_o, outputs = ParallelWrapper.assemble() 
 
  
 
class SLABPLOT(LINEARPLOT):
    '''
    Plot slab morphology
    Inputs:
        -
    Returns:
        -
    '''
    def __init__(self, _name):
        LINEARPLOT.__init__(self, _name)
 
 
    def PlotTrenchPosition(self, case_dir, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        time_interval_for_slab_morphology = kwargs.get("time_interval", 1.0e6)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
        trench_edge_y = float(Visit_Options.options['TRENCH_EDGE_Y_FULL'])
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        assert(len(available_pvtu_snapshots) > 1) # assert this is an array
        n_snapshots = len(available_pvtu_snapshots)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # output the availabe pvtu snapshots
        
        # make directory to save images 
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_img_dir = os.path.join(img_dir, "morphology")
        if not os.path.isdir(morph_img_dir):
            os.mkdir(morph_img_dir)

        # 1. trench position
        n_in_group = 5
        n_column = 3
        n_group = int(np.ceil(1.0 * n_snapshots / n_in_group))
        n_raw = int(np.ceil(1.0 * n_group / n_column))

        fig = plt.figure(tight_layout=True, figsize=(5 * n_column, 5 * n_raw)) 
        gs = gridspec.GridSpec(n_raw, n_column) 
        for i in range(n_group):
            col = i % n_column
            raw = i // n_column
            ax = fig.add_subplot(gs[raw, col]) 
            for j in range(n_in_group+1):
                # j value is chosen to connect the last one in the
                # previous group and the first one in the current group
                i_in_snapshots = i * n_in_group + j
                if i_in_snapshots < n_snapshots:
                    vtu_snapshot = available_pvtu_snapshots[i_in_snapshots]
                    _time = time_interval_for_slab_morphology * i_in_snapshots
                    normalizer = [ float(j)/(n_in_group) for i in range(n_in_group+1)] 
                    colors = cm.rainbow(normalizer) 
                    _color = colors[j]
                    ax = self.PlotTrenchPositionStep(case_dir, vtu_snapshot, axis=ax, color=_color, label="Time = %.2f Myr" % (_time / 1e6))
            ax.legend()
            ax.set_xlim([trench_initial_position/1e3 - 1000, trench_initial_position/1e3 + 1000]) # set the x axis to center around the inital trench
        # save figure
        fileout = os.path.join(morph_img_dir, "trench_history.png")
        fig.savefig(fileout)
        print("image output: ", fileout)

        # 2. trench velocity
        fig, ax= plt.subplots(tight_layout=True, figsize=(8, 5)) 
        # use the steps between the 0th and 1th snapshot as the interval for calculating the velocities
        step_for_derivatives = available_pvtu_snapshots[1] - available_pvtu_snapshots[0]
        y_queries_for_trench_velocities = np.arange(0.0, trench_edge_y - 100e3, 200e3)
        silence = False
        for i in range(y_queries_for_trench_velocities.size):
            if i > 0:
                silence = True
            y_query = y_queries_for_trench_velocities[i]
            ax = self.PlotTrenchVelocity(case_dir, time_interval_for_slab_morphology,\
                                        axis=ax, step_for_derivatives=step_for_derivatives, y_query=y_query, silence=silence)
        ax.legend()
        # save image
        fileout = os.path.join(morph_img_dir, "trench_velocities.png")
        fig.savefig(fileout)
        print("image output: ", fileout)

    def PlotTrenchPositionAnimation(self, case_dir, vtu_step_list, this_vtu_step, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        time_interval_for_slab_morphology = kwargs.get("time_interval", 1.0e6)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
        trench_edge_y = float(Visit_Options.options['TRENCH_EDGE_Y_FULL'])
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        assert(len(available_pvtu_snapshots) > 1) # assert this is an array
        n_snapshots = len(available_pvtu_snapshots)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # output the availabe pvtu snapshots
        
        # make directory to save images 
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_img_dir = os.path.join(img_dir, "morphology")
        if not os.path.isdir(morph_img_dir):
            os.mkdir(morph_img_dir)
        
        # 1. trench position
        n_episode = len(vtu_step_list)
        fig, ax = plt.subplots(tight_layout=True, figsize=(5, 5))
        normalizer = [float(i)/(n_episode+1) for i in range(n_episode)] 
        colors = cm.rainbow(normalizer) 
        # 
        for i in range(n_episode):
            vtu_snapshot = vtu_step_list[i] + int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT'])
            _time, _ = Visit_Options.get_time_and_step(vtu_step_list[i])
            _color = colors[i]
            ax = self.PlotTrenchPositionStep(case_dir, vtu_snapshot, axis=ax, color=_color, label="Time = %.2f Myr" % (_time / 1e6))
            ax.legend()
            ax.set_xlim([trench_initial_position/1e3 - 1000, trench_initial_position/1e3 + 1000]) # set the x axis to center around the inital trench
        #
        this_vtu_snapshot = this_vtu_step + int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT'])
        this_time, _ = Visit_Options.get_time_and_step(this_vtu_step)
        ax = self.PlotTrenchPositionStep(case_dir, this_vtu_snapshot, axis=ax, color=_color, label="Time = %.2f Myr" % (this_time / 1e6), line_type="-")
        fileout = os.path.join(morph_img_dir, "trench_history_animation_t%.4e.png" % this_time)
        fig.savefig(fileout)
        print("image output: ", fileout)

    
    def PlotTrenchPositionEpisodes(self, case_dir, episodes, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results with episodes
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        time_interval_for_slab_morphology = kwargs.get("time_interval", 1.0e6)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
        trench_edge_y = float(Visit_Options.options['TRENCH_EDGE_Y_FULL'])
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        assert(len(available_pvtu_snapshots) > 1) # assert this is an array
        n_snapshots = len(available_pvtu_snapshots)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # output the availabe pvtu snapshots
        
        # make directory to save images 
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_img_dir = os.path.join(img_dir, "morphology")
        if not os.path.isdir(morph_img_dir):
            os.mkdir(morph_img_dir)

        # 1. trench position
        n_column = 2
        n_episodes = len(episodes)
        n_raw = int(np.ceil(1.0 * n_episodes / n_column))

        fig = plt.figure(tight_layout=True, figsize=(5 * n_column, 5 * n_raw)) 
        gs = gridspec.GridSpec(n_raw, n_column) 
        for i in range(n_episodes):
            col = i % n_column
            raw = i // n_column
            ax = fig.add_subplot(gs[raw, col]) 
            episode = episodes[i]
            assert(len(episode) == 2)
            for i_in_snapshots in range(episode[0], episode[1]):
                if i_in_snapshots < n_snapshots:
                    vtu_snapshot = available_pvtu_snapshots[i_in_snapshots]
                    _time = time_interval_for_slab_morphology * i_in_snapshots
                    normalizer = [ float(i_in_snapshots - episode[0])/(episode[1] - episode[0] -1) for i in range(episode[1] - episode[0])] 
                    colors = cm.rainbow(normalizer) 
                    _color = colors[i_in_snapshots - episode[0]]
                    ax = self.PlotTrenchPositionStep(case_dir, vtu_snapshot, axis=ax, color=_color, label="Time = %.2f Myr" % (_time / 1e6))
            ax.legend()
            ax.set_xlim([trench_initial_position/1e3 - 1000, trench_initial_position/1e3 + 1000]) # set the x axis to center around the inital trench
        # save figure
        fileout = os.path.join(morph_img_dir, "trench_history_episodes.png")
        fig.savefig(fileout)
        print("image output: ", fileout)


    def PlotTrenchPositionStep(self, case_dir, vtu_snapshot, **kwargs):
        '''
        plot trench position for a single step
        '''
        _color = kwargs.get("color", None)
        _label = kwargs.get("label", None)
        line_type = kwargs.get("line_type", '.')
        filein = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % vtu_snapshot)
        Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found." % filein)
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        data = np.loadtxt(filein)
        xs = data[:, 0]
        ys = data[:, 1]
        ax.plot(xs/1e3, ys/1e3, line_type, color=_color, label=_label)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        return ax
    
    
    def PlotTrenchVelocity(self, case_dir, time_interval_for_slab_morphology, **kwargs):
        '''
        plot trench position for a single step
        '''
        _color = kwargs.get("color", None)
        _label = kwargs.get("label", None)
        silence = kwargs.get("silence", False)
        step_for_derivatives = kwargs.get('step_for_derivatives', 5)
        y_query = kwargs.get("y_query", 0.0)

        # initiation
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        n_snapshots = len(available_pvtu_snapshots)
        # initiate plot
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
       
        # derive the trench velocities
        V_trs = []
        ts = []
        for n in range(n_snapshots-1):
            vtu_snapshot = available_pvtu_snapshots[n]
            filein0 = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % vtu_snapshot)
            filein1 = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % (vtu_snapshot+step_for_derivatives))
            if not os.path.isfile(filein0):
                if not silence:
                    warnings.warn('PlotTrenchVelocity: File %s is not found' % filein0, UserWarning)
                continue
            elif not os.path.isfile(filein1):
                if not silence:
                    warnings.warn('PlotTrenchVelocity: File %s is not found' % filein1, UserWarning)
                continue
            # Utilities.my_assert(os.path.isfile(filein0), FileNotFoundError, "%s doens't exist" % filein0)
            # Utilities.my_assert(os.path.isfile(filein1), FileNotFoundError, "%s doens't exist" % filein1)
            _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
            _time1, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot+step_for_derivatives)
            dt = _time1 - _time
            data0 = np.loadtxt(filein0)
            xs0 = data0[:, 0]
            ys0 = data0[:, 1]
            x_tr0 = np.interp(y_query, ys0, xs0)
            data1 = np.loadtxt(filein1)
            xs1 = data1[:, 0]
            ys1 = data1[:, 1]
            x_tr1 = np.interp(y_query, ys1, xs1)
            V_tr0 = (x_tr1 - x_tr0) / dt
            V_trs.append(V_tr0)
            ts.append((_time1 + _time) / 2.0)
        V_trs = np.array(V_trs)
        ts = np.array(ts)
        ax.plot(ts / 1e6, V_trs * 100.0, label="y = %.2f km" % (y_query / 1e3))
        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Velocity (cm/yr)')
        ax.grid()
        return ax


class PLOT_COMBINE_SLAB_MORPH(PLOT_COMBINE):
    '''
    Combine results from slab morphology
    '''
    def __init__(self, case_paths):
        PLOT_COMBINE.__init__(self, case_paths)
        UnitConvert = Utilities.UNITCONVERT()
        self.MorphPlotter = SLABPLOT("plot_slab")
        pass

    def __call__(self, width, output_dir, time_range, tp_range, sd_range, **kwargs):
        '''
        perform combination
        Inputs:
            sizes: (list of 2) - size of the plot
            output_dir: directory to output to
        '''
        _name = "combine_morphology"
        _title = "Comparing slab morphology results"
        color_method = kwargs.get('color_method', 'generated')
        dump_color_to_json = kwargs.get('dump_color_to_json', None)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

def PlotTrenchDifferences(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    _color = kwargs.get("color", None)
    _label = kwargs.get("label", None)
    silence = kwargs.get("silence", False)
    y_query = kwargs.get("y_query", 0.0)
    ax_twinx = kwargs.get("axis_twinx", None)

    # initiation
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)
    # initiate plot
    ax = kwargs.get('axis', None)
    if ax == None:
        raise ValueError("Not implemented")
   
    # derive the trench locations and interpolate to the query points
    # derive the depth of the slab tip
    dx_trs = []
    depths = []
    ts = []
    for n in range(n_snapshots-1):
        vtu_snapshot = available_pvtu_snapshots[n]
        filein0 = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % vtu_snapshot)
        filein1 = os.path.join(case_dir, "vtk_outputs", "center_profile_%05d.txt" % vtu_snapshot)
        if not os.path.isfile(filein0):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein0, UserWarning)
            continue
        if not os.path.isfile(filein1):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein1, UserWarning)
            continue
        _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
        # get the trench locations
        data0 = np.loadtxt(filein0)
        xs0 = data0[:, 0]
        ys0 = data0[:, 1]
        x_tr0 = np.interp(y_query, ys0, xs0)
        dx_trs.append(x_tr0)
        # get the slab tip depth
        data1 = np.loadtxt(filein1)
        zs1 = data1[:, 2]
        depth = float(Visit_Options.options["BOX_THICKNESS"]) - np.min(zs1)
        depths.append(depth)
        ts.append(_time)
    dx_trs = np.array(dx_trs)
    depths = np.array(depths)
    ts = np.array(ts)
    # plot the differences by subtracting the initial value
    ax.plot(ts / 1e6, (dx_trs - dx_trs[0]) / 1e3, label="y = %.2f km" % (y_query / 1e3), color=_color) # trench position
    if y_query < 1e-6:
        # only plot the depth at the center
        ax_twinx.plot(ts / 1e6, depths / 1e3, "--", color=_color) # depth
    return ax


def PlotSlabDipAngle(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    _color = kwargs.get("color", None)
    _label = kwargs.get("label", None)
    silence = kwargs.get("silence", False)
    y_query = kwargs.get("y_query", 0.0)
    ax_twinx = kwargs.get("axis_twinx", None)
    pin_depth = kwargs.get("pin_depth", 100e3)
    crust_only = kwargs.get("crust_only", 0)

    # initiation
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)
    # initiate plot
    ax = kwargs.get('axis', None)
    if ax == None:
        raise ValueError("Not implemented")
   
    # derive the trench locations and interpolate to the query points
    # derive the depth of the slab tip
    dips = []
    depths = []
    ts = []
    for n in range(n_snapshots-1):
        vtu_snapshot = available_pvtu_snapshots[n]
        # mark the composition used for the envelop in the filename
        if crust_only == 1:
            filename0 = "trench_b_%05d.txt" % vtu_snapshot
            filename1 = "center_profile_%05d.txt" % vtu_snapshot
            filename2 = "trench_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
        elif crust_only == 2:
            filename0 = "trench_l_b_%05d.txt" % vtu_snapshot
            filename1 = "center_profile_l_%05d.txt" % vtu_snapshot
            filename2 = "trench_l_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
        else:
            filename0 = "trench_lu_b_%05d.txt" % vtu_snapshot
            filename1 = "center_profile_lu_%05d.txt" % vtu_snapshot
            filename2 = "trench_lu_b_d%.2fkm_%05d.txt" % (pin_depth/1e3, vtu_snapshot)
        filein0 = os.path.join(case_dir, "vtk_outputs", filename0)
        filein1 = os.path.join(case_dir, "vtk_outputs", filename1)
        filein2 = os.path.join(case_dir, "vtk_outputs", filename2)
        if not os.path.isfile(filein0):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein0, UserWarning)
            continue
        if not os.path.isfile(filein1):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein1, UserWarning)
            continue
        _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
        # get the trench locations
        data0 = np.loadtxt(filein0)
        xs0 = data0[:, 0]
        ys0 = data0[:, 1]
        zs0 = data0[:, 2]
        x_tr0 = np.interp(y_query, ys0, xs0)
        z_tr0 = np.interp(y_query, ys0, zs0)
        # get the slab tip depth
        data1 = np.loadtxt(filein1)
        zs1 = data1[:, 2]
        depth = float(Visit_Options.options["BOX_THICKNESS"]) - np.min(zs1)
        depths.append(depth)
        ts.append(_time)
        # get the curve at a certain depth
        data2 = np.loadtxt(filein2)
        xs2 = data2[:, 0]
        ys2 = data2[:, 1]
        zs2 = data2[:, 2]
        x_tr2 = np.interp(y_query, ys2, xs2)
        z_tr2 = np.interp(y_query, ys2, zs2)
        dip = - math.atan((z_tr0 - z_tr2)/(x_tr0 - x_tr2))
        print("vtu_snapshot = %d, x_tr0 = %.4e, z_tr0 = %.4e, x_tr2 = %.4e, z_tr2 = %.4e, dip = %.4e" % (vtu_snapshot, x_tr0, z_tr0, x_tr2, z_tr2, dip))  # debug
        dips.append(dip)
    dips = np.array(dips)
    depths = np.array(depths)
    ts = np.array(ts)
    # apply a spline
    spline = UnivariateSpline(ts/1e6, dips, s=0) # s=0 means interpolation
    ts_new = np.linspace(ts.min(), ts.max(), 1000)
    dips_splined = spline(ts_new/1e6)
    # plot the differences by subtracting the initial value
    ax.plot(ts / 1e6, dips * 180.0 / np.pi, label="y = %.2f km" % (y_query / 1e3), color=_color) # trench position
    ax.plot(ts_new / 1e6, dips_splined * 180.0 / np.pi, "--", label="y = %.2f km (splined)" % (y_query / 1e3), color=_color) # trench position
    if y_query < 1e-6:
        # only plot the depth at the center
        ax_twinx.plot(ts / 1e6, depths / 1e3, "--", color=_color) # depth
    return ax


def PlotSlabDepth(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    _color = kwargs.get("color", None)
    _label = kwargs.get("label", None)
    silence = kwargs.get("silence", False)
    y_query = kwargs.get("y_query", 0.0)

    # initiation
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)
    # initiate plot
    ax = kwargs.get('axis', None)
    if ax == None:
        raise ValueError("Not implemented")
   
    # derive the depth of the slab tip
    depths = []
    ts = []
    for n in range(n_snapshots-1):
        vtu_snapshot = available_pvtu_snapshots[n]
        filein1 = os.path.join(case_dir, "vtk_outputs", "center_profile_%05d.txt" % vtu_snapshot)
        if not os.path.isfile(filein1):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein1, UserWarning)
            continue
        _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
        # get the slab tip depth
        data1 = np.loadtxt(filein1)
        zs1 = data1[:, 2]
        depth = float(Visit_Options.options["BOX_THICKNESS"]) - np.min(zs1)
        depths.append(depth)
        ts.append(_time)
    depths = np.array(depths)
    ts = np.array(ts)
    if y_query < 1e-6:
        # only plot the depth at the center
        ax.plot(ts / 1e6, depths / 1e3, "--", color=_color) # depth
    return ax


def CenterProfileAnalyze(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    silence = kwargs.get("silence", 0)
    # get the index in an array by a given value
    IndexByValue = lambda array_1d, val: np.argmin(abs(array_1d - val))
    # initiation
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)

    # derive the depth of the slab tip
    depths = []
    ts = []
    steps = []
    for n in range(n_snapshots-1):
        vtu_snapshot = available_pvtu_snapshots[n]
        filein1 = os.path.join(case_dir, "vtk_outputs", "center_profile_%05d.txt" % vtu_snapshot)
        if not os.path.isfile(filein1):
            if not silence:
                warnings.warn('PlotTrenchVelocity: File %s is not found' % filein1, UserWarning)
            continue
        _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
        # get the slab tip depth
        data1 = np.loadtxt(filein1)
        zs1 = data1[:, 2]
        depth = float(Visit_Options.options["BOX_THICKNESS"]) - np.min(zs1)
        depths.append(depth)
        ts.append(_time)
        steps.append(step)
    depths = np.array(depths)
    ts = np.array(ts)

    # time of slab tip reaching 660 km and the index in the list
    sfunc = interp1d(depths, ts, assume_sorted=True)
    t660 = sfunc(660e3)
    i660 = IndexByValue(ts, t660)
    step660 = steps[i660]

    return t660, step660


# functions
def GetVtkCells2d(n_x, n_z):
    '''
    Get a vtk cells from number of points along x, z in 2 dimensions
    '''
  
    cells = vtk.vtkCellArray()
  
    for ix in range(n_x-1):
        for jz in range(n_z-1):
          cells.InsertNextCell(4)
          # cell = vtk.vtkIdType()
          cells.InsertCellPoint(ix*n_z + jz+1)
          cells.InsertCellPoint(ix*n_z + jz)
          cells.InsertCellPoint((ix+1)*n_z + jz)
          cells.InsertCellPoint((ix+1)*n_z + jz+1)
  
    return cells


def Interpolate3dVtkCase(case_dir, vtu_snapshot, **kwargs):
    '''
    interpolate a 3d vtk output from a case
    Inputs:
        kwargs:
            interval - this determines the interval of the slices
    '''
    interval = kwargs.get("interval", 10e3)
    file_extension = kwargs.get("file_extension", "vtp")
    n0 = kwargs.get("n0", 800)
    n1 = kwargs.get("n1", 100)

    assert(os.path.isdir(case_dir))

    ha_file = os.path.join(case_dir, "output", "depth_average.txt")
    assert(os.path.isfile(ha_file))
    
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Ri = Visit_Options.options['INNER_RADIUS']
    Xmax = Visit_Options.options['XMAX']
    print("geometry: %s, Ro/Depth: %f, Xmax: %f" % (geometry, Ro, Xmax))
    
    # load the class for processing the depth average output
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.Import(ha_file)

    # parse file location, time and step
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found" % filein)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)

    # parse from depth average files 
    Utilities.my_assert(_time != None, ValueError, "\"time\" is a requried input if \"ha_file\" is presented")
    Tref_func = DepthAverage.GetInterpolateFunc(_time, "temperature")
    density_ref_func = DepthAverage.GetInterpolateFunc(_time, "adiabatic_density")

    # initiate vtk readers
    reader = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(filein)

    # read data
    start = time.time() # time this operation
    reader.Update()
    grid = reader.GetOutput()
    data_set = reader.GetOutputAsDataSet()
    points = grid.GetPoints()
    cells = grid.GetCells()
    point_data = data_set.GetPointData()
    end = time.time()
    print("Read data takes %.2f s" % (end-start))

    # get points and point data 
    points_np = vtk_to_numpy(points.GetData())
    T_np = vtk_to_numpy(point_data.GetArray("T"))

    # new mesh
    N = n0 * n1
    target_points_np = np.zeros((N, 3))
    target_cells_vtk = GetVtkCells2d(n0, n1)
    mask=None
    if geometry == "box":
        y0 = 1.0
        for i0 in range(n0):
            for j1 in range(n1):
                ii = i0 * n1 + j1
                target_points_np[ii, 0] = Xmax * i0 / (n0 - 1)
                target_points_np[ii, 1] = y0
                target_points_np[ii, 2] = Ro * j1 / (n1 - 1) # Ro and depth are the same in this case
        mask = (points_np[:, 1] > (y0 - interval)) & (points_np[:, 1] < (y0 + interval))
    elif geometry == "chunk":
        z0 = 1.0
        for i0 in range(n0):
            for j1 in range(n1):
                # note we use theta = 0.0 here, but z0 = small value, this is to ensure a successful slice
                # of the chunk geometry
                ii = i0 * n1 + j1
                phi = i0 / n0 * Xmax * np.pi / 180.0
                r = Ro * (j1 - 0.0) / (n1 - 0.0) + Ri * (j1 - n1)/ (0.0 - n1)  
                slice_x, slice_y, _  = Utilities.ggr2cart(0.0, phi, r) 
                target_points_np[ii, 0] = slice_x
                target_points_np[ii, 1] = slice_y
                target_points_np[ii, 2] = z0
        mask = (points_np[:, 2] > (z0 - interval)) & (points_np[:, 2] < (z0 + interval))

    # apply a mask before interpolate
    print("Interval = %.2e, Adjacent points number = %d" % (interval, points_np[mask].shape[0]))

    # Resample the data from the source mesh to the target mesh
    resampled_data = griddata(points_np[mask], T_np[mask], target_points_np, method='linear')

    # path to output
    odir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(odir):
        os.mkdir(odir)
    fileout=os.path.join(odir, "center_slice-%05d.%s" % (vtu_snapshot, file_extension))

    # output to vtp file
    o_poly_data = vtk.vtkPolyData()  # initiate poly daa
    # insert points
    target_points_vtk = vtk.vtkPoints()
    for i in range(target_points_np.shape[0]):
      target_points_vtk.InsertNextPoint(target_points_np[i, 0], target_points_np[i, 1], target_points_np[i, 2])
    o_poly_data.SetPoints(target_points_vtk) # insert points
    o_poly_data.SetPolys(target_cells_vtk)
    # insert data
    resampled_T_vtk = numpy_to_vtk(resampled_data, deep=1)
    resampled_T_vtk.SetName('T')
    o_poly_data.GetPointData().SetScalars(resampled_T_vtk)

    # initiate writer and output 
    if file_extension == 'txt':
        if resampled_data.ndim == 1:
            resampled_data = resampled_data.reshape(resampled_data.size, 1)
        odata = np.concatenate((target_points_np, resampled_data), axis=1)
        np.savetxt(fileout, odata)
        # writer = vtk.vtkSimplePointsWriter()
    elif file_extension == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileout)
        writer.SetInputData(o_poly_data)
        # writer.SetFileTypeToBinary()  # try this later to see if this works
        writer.Update()
        writer.Write()
    else:
        raise NotImplementedError()
    
    print("Saved file %s" % fileout)


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
    parser.add_argument('-d', '--depth', type=float,
                        default=150e3,
                        help='Depth of the cross section')
    parser.add_argument('-vs', '--vtu_step', type=int,
                        default=0,
                        help='vtu_step')
    parser.add_argument('-vss', '--vtu_snapshot', type=int,
                        default=0,
                        help='vtu_snapshot')
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=1.0e6,
                        help='Time interval, affecting the time steps to visualize')
    parser.add_argument('-seix', '--slab_envelop_interval_y', type=float,
                        default=20e6,
                        help='Interval along y axis to sort out the trench locations')
    parser.add_argument('-seiz', '--slab_envelop_interval_z', type=float,
                        default=20e6,
                        help='Interval along z axis to sort out the trench locations')
    parser.add_argument('-ssc', '--slab_shallow_cutoff', type=float,
                        default=40e6,
                        help='Minimum depth along z axis to sort out the trench locations')
    parser.add_argument('-co', '--crust_only', type=int,
                        default=0,
                        help='If we only use the crustal composition to sort out the trench locations')
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
        # include_force_balance = True
        include_force_balance = False
        SlabMorphology(arg.inputs, int(arg.vtu_snapshot), rewrite=1,\
            slab_envelop_interval_y=arg.slab_envelop_interval_y, slab_envelop_interval_z=arg.slab_envelop_interval_z,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only, include_force_balance=include_force_balance)
    elif _commend == 'cross_section_at_depth':
        SlabMorphology(arg.inputs, int(arg.vtu_snapshot), rewrite=1,\
            slab_envelop_interval_y=arg.slab_envelop_interval_y, slab_envelop_interval_z=arg.slab_envelop_interval_z,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only, horizontal_velocity_depths=[arg.depth])
    elif _commend == 'morph_case':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=0, time_interval=arg.time_interval,\
            slab_envelop_interval_y=arg.slab_envelop_interval_y, slab_envelop_interval_z=arg.slab_envelop_interval_z,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only)
    elif _commend == 'morph_case_parallel':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=1, time_interval=arg.time_interval,\
            slab_envelop_interval_y=arg.slab_envelop_interval_y, slab_envelop_interval_z=arg.slab_envelop_interval_z,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only, use_parallel=True)
    elif _commend == 'plot_trench':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotTrenchPosition(arg.inputs, time_interval=arg.time_interval)
    elif _commend == 'plot_trench_by_episodes':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        episodes = [[0, 5], [5, 20], [20, 30]]
        SlabPlot.PlotTrenchPositionEpisodes(arg.inputs, episodes, time_interval=arg.time_interval)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()