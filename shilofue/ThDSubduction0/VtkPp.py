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
from posixpath import split
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
from shilofue.TwoDSubduction0.VtkPp import get_theta, get_dip
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
        self.slab_envelop_interval_d = kwargs.get("slab_envelop_interval_d", 10e3)
        self.slab_envelop_interval_w = kwargs.get("slab_envelop_interval_w", 10e3)
        self.slab_envelop_point_list0 = []
        self.slab_envelop_point_list1 = []
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
        start = time.time()
        print("PrepareSlabByPoints")
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
        end = time.time()
        print("\tInitiation, takes %.2f s" % (end-start))
        start = end
        # add cells by composition
        # todo_3d_chunk
        min_r = self.Ro
        slab_w_max = 0.0
        for i in range(self.i_poly_data.GetNumberOfPoints()):
            x = points[i][0]
            y = points[i][1]
            z = points[i][2]
            r, w, _ = get_slab_dimensions_3(x, y, z, self.Ro, self.is_chunk)
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
                if w > slab_w_max:
                    slab_w_max = w
        self.slab_depth = self.Ro - min_r  # cart
        end = time.time()
        print("\tLook for slab points, takes %.2f s" % (end-start))
        print("\t%d points found in the subducting slab" % len(self.slab_points))
        start = end
        # prepare the slab internal points
        # The length, width and depth are used to determine 
        # the location of points within the slab consistently
        # between the cartesian and the chunk geometry
        total_en_interval_d = int((self.slab_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval_d + 1)
        total_en_interval_w = int((slab_w_max) // self.slab_envelop_interval_w + 1)
        id_en_pin_depth = int((pin_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval_d)
        slab_en_point_lists = [ [] for i in range(total_en_interval_d * total_en_interval_w) ]
        for id in self.slab_points:
            x = points[id][0] # first, separate cells into intervals
            y = points[id][1]
            z = points[id][2]
            r, w, _ = get_slab_dimensions_3(x, y, z, self.Ro, self.is_chunk)
            # todo_3d_chunk
            id_en_d =  int(np.floor(
                                  (self.Ro - r - self.slab_shallow_cutoff)/
                                  self.slab_envelop_interval_d))# id in the envelop list
            id_en_w =  int(np.floor(w / self.slab_envelop_interval_w))
            id_en = id_en_w * total_en_interval_d + id_en_d
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
        end = time.time()
        print("\tCategorize slab points, takes %.2f s" % (end-start))
        start = end
        # extract the slab envelop
        # The internal points are processed by their length dimension
        # consistently between the cartesian and the chunk geometry.
        # The slab envelops have the min and max value in the length dimension
        for id_en in range(len(slab_en_point_lists)):
            l_min = 0.0  # then, loop again by intervals to look for a
            l_max = 0.0  # max w and a min w for each interval
            point_list = slab_en_point_lists[id_en]
            if len(point_list) == 0:
                continue  # make sure we have some point
            is_first = True
            id_min = -1
            id_max = -1
            for id in point_list:
                x = points[id][0]
                y = points[id][1]
                z = points[id][2]
                _, _, l = get_slab_dimensions_3(x, y, z, self.Ro, self.is_chunk)
                # todo_3d_chunk
                if is_first:
                    id_min = id
                    id_max = id
                    l_min = l
                    l_max = l
                    is_first = False
                else:
                    if l < l_min:
                        id_min = id
                        l_min = l
                    if l > l_max:
                        id_max = id
                        l_max = l
            ''' DEBUG OUTPUTS - output the points on the top and bottom of the envelops
            id_w = int(id_en // total_en_interval_d)
            id_en_d = int(id_en % total_en_interval_d)
            if id_w == 0: 
                x_min, y_min, z_min = points[id_min][0], points[id_min][1], points[id_min][2]
                r_min, th_min, ph_min = Utilities.cart2sph(x_min,y_min,z_min)
                x_max, y_max, z_max = points[id_max][0], points[id_max][1], points[id_max][2]
                r_max, th_max, ph_max = Utilities.cart2sph(x_max,y_max,z_max)
                print("%d: (%.4e, %.4e, %.4e), (%.4e, %.4e, %.4e)" % (id_en_d, x_min, y_min, z_min, r_min, th_min, ph_min))
                print("%d: (%.4e, %.4e, %.4e), (%.4e, %.4e, %.4e)" % (id_en_d, x_max, y_max, z_max, r_max, th_max, ph_max))
            '''
            self.slab_envelop_point_list0.append(id_min)  # first half of the envelop
            self.slab_envelop_point_list1.append(id_max)  # second half of the envelop
            if id_en < total_en_interval_d:
                # record center profile
                # todo_center
                if id_max > 0:
                    center_profile_x.append(points[id_max][0])
                    center_profile_y.append(points[id_max][1])
                    center_profile_z.append(points[id_max][2])
            if id_en % total_en_interval_d == 0:
                # record trench coordinates
                id_en_w = id_en // total_en_interval_d
                if id_max > 0:
                    trench_coords_x.append(points[id_max][0])
                    trench_coords_y.append(points[id_max][1])
                    trench_coords_z.append(points[id_max][2])
                if id_min > 0:
                    trench_b_coords_x.append(points[id_min][0])
                    trench_b_coords_y.append(points[id_min][1])
                    trench_b_coords_z.append(points[id_min][2])
            if get_pin_depth:
                if id_en % total_en_interval_d == id_en_pin_depth:
                    # get the pinned x, y and z at a pinned depth
                    if id_max > 0:
                        pinD_coords_x.append(points[id_max][0])
                        pinD_coords_y.append(points[id_max][1])
                        pinD_coords_z.append(points[id_max][2])
                    if id_min > 0:
                        pinD_b_coords_x.append(points[id_min][0])
                        pinD_b_coords_y.append(points[id_min][1])
                        pinD_b_coords_z.append(points[id_min][2])
        end = time.time()
        print("\tExtract slab envelops (top and bottom), takes %.2f s" % (end-start))
        start = end
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
        end = time.time()
        print("\tWrap up, takes %.2f s" % (end-start))
        start = end

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
            r = VtkPp.get_r3(x, y, z, self.is_chunk)
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
    Wrapper for using the PVTK class to analyze one step.

    Inputs:
        case_dir (str): Path to the case directory.
        vtu_snapshot (int): Step number for VTU outputs.
        kwargs (dict): Additional parameters for the analysis:
            - slab_envelop_interval_w (float, default=10e3): Interval along the y-axis for sorting trench locations.
            - slab_envelop_interval_d (float, default=10e3): Interval along the z-axis for sorting trench locations.
            - slab_shallow_cutoff (float, default=70e3): Minimum depth along the z-axis for sorting trench locations.
            - include_force_balance (bool, default=False): Whether to include force balance calculations.
            - crust_only (int, default=0): Specifies which crustal compositions to consider for trench sorting.
            - horizontal_velocity_depths (list, default=[]): Depths at which to investigate the velocity field.
            - pin_depth (float, default=100e3): Depth at which points are pinned.
            - output (str, optional): Output path for VTK files.
    
    Description:
        - Asserts that the input solution file exists in the specified case directory.
        - Sets or retrieves default parameters for trench location analysis.
        - Initializes the VISIT_OPTIONS class to interpret geometry and other setup information.
        - Initializes the VTKP object and prepares data for analysis.
        - Constructs polydata and processes slab composition based on specified crust types.
        - Exports data, including slab internal points, envelopes, trench positions, center profile, and pinned depth data.
        - If specified, includes force balance calculations and exports buoyancy and density arrays.
        - Outputs trench coordinates, slab envelopes, and other related data to specified files.
    
    Returns:
        tuple: The vtu_step used and an empty string (for compatibility).
    '''
    
    # Construct the file path for the input solution file and assert it exists
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found" % filein)
    
    # Retrieve or set default values for various parameters
    slab_envelop_interval_w = kwargs.get("slab_envelop_interval_w", 10e3)
    slab_envelop_interval_d = kwargs.get("slab_envelop_interval_d", 10e3)
    slab_shallow_cutoff = kwargs.get("slab_shallow_cutoff", 70e3)
    include_force_balance = kwargs.get("include_force_balance", False)
    crust_only = kwargs.get("crust_only", 0)
    horizontal_velocity_depths = kwargs.get("horizontal_velocity_depths", [])
    pin_depth = kwargs.get("pin_depth", 100e3)
    appendix = "_%05d" % vtu_snapshot
    
    # Set the output path and create the directory if it does not exist
    output_path = kwargs.get("output", os.path.join(case_dir, 'vtk_outputs'))
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    
    # Initialize VISIT_OPTIONS and interpret the case's geometry
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro = Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    print("geometry: %s, Ro: %f, Xmax: %f" % (geometry, Ro, Xmax))
    
    # Initialize the VTKP object with appropriate parameters
    if include_force_balance:
        ha_file = os.path.join(case_dir, "output", "depth_average.txt")
        assert(os.path.isfile(ha_file))
    else:
        ha_file = None
    
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax, slab_shallow_cutoff=slab_shallow_cutoff,
                slab_envelop_interval_w=slab_envelop_interval_w, slab_envelop_interval_d=slab_envelop_interval_d,
                ha_file=ha_file, time=_time)
    VtkP.ReadFile(filein)
    
    # Construct polydata and record the processing time
    field_names = ['T', 'p', 'density', 'sp_upper', 'sp_lower']
    start = time.time()
    VtkP.ConstructPolyData(field_names, include_cell_center=False)
    end = time.time()
    print("Construct polydata, takes %.2f s" % (end - start))
    start = end
    
    # Prepare slab composition data points based on crust_only flag
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
    
    # Export slab internal points and handle force balance if needed
    slab_point_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data, VtkP.slab_points)
    if include_force_balance:
        # Mark the composition used for the envelope in the filename
        if crust_only == 1:
            filename = 'slab.vtp'
        elif crust_only == 2:
            filename = 'slab_l.vtp'
        else:
            filename = 'slab_lu.vtp'
        fileout = os.path.join(output_path, filename)
        # Compute buoyancy and prepare output data
        VtkP.ComputeBuoyancy()
        buoyancy_vtk_array = numpy_to_vtk(VtkP.buoyancies)
        buoyancy_vtk_array.SetName("buoyancy")
        o_poly_data = vtk.vtkPolyData()
        o_poly_data.SetPoints(slab_point_grid.GetPoints())
        o_poly_data.GetPointData().SetScalars(buoyancy_vtk_array)
        # Add density arrays and export data
        densities = vtk_to_numpy(VtkP.i_poly_data.GetPointData().GetArray('density'))
        densities_slab = np.zeros(len(VtkP.slab_points))
        for i in range(len(VtkP.slab_points)):
            id_en = VtkP.slab_points[i]
            densities_slab[i] = densities[id_en]
        density_vtk_array = numpy_to_vtk(densities_slab)
        density_vtk_array.SetName("density")
        o_poly_data.GetPointData().AddArray(density_vtk_array)
        densities_ref = np.zeros(len(VtkP.slab_points))
        for i in range(len(VtkP.slab_points)):
            id_en = VtkP.slab_points[i]
            x, y, z = VtkP.i_poly_data.GetPoints().GetData()[id_en]
            r = VtkPp.get_r3(x, y, z, VtkP.is_chunk)
            depth = Ro - r
            densities_ref[i] = VtkP.density_ref_func(depth)
        densities_ref_vtk_array = numpy_to_vtk(densities_ref)
        densities_ref_vtk_array.SetName("density_ref")
        o_poly_data.GetPointData().AddArray(densities_ref_vtk_array)
        VtkPp.ExportPolyData(o_poly_data, fileout)
    else:
        # Export data without force balance calculations
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
        # Mark the composition used for the envelope in the filename
        filename = "slab_env0" + appendix + ".vtu"
    elif crust_only == 2:
        filename = "slab_env0_l" + appendix + ".vtu"
    else:
        filename = "slab_env0_lu" + appendix + ".vtu"
    fileout = os.path.join(output_path, filename)
    slab_env0_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data, VtkP.slab_envelop_point_list0)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_env0_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab envelope0 (smaller x) points: %s" % fileout)

    # 2b Top envelope
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
        filename = "slab_env1" + appendix + ".vtu"
    elif crust_only == 2:
        filename = "slab_env1_l" + appendix + ".vtu"
    else:
        filename = "slab_env1_lu" + appendix + ".vtu"
    fileout = os.path.join(output_path, filename)
    slab_env1_grid = VtkPp.ExportPointGridFromPolyData(VtkP.i_poly_data, VtkP.slab_envelop_point_list1)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_env1_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    print("File output for slab envelope1 (bigger x) points: %s" % fileout)
    end = time.time()
    print("Write slab points, takes %.2f s" % (end - start))
    start = end

    # Prepare outputs for trench positions
    # a. Top trench positions
    trench_coords_x, trench_coords_y, trench_coords_z = VtkP.ExportSlabInfo()
    outputs = "# trench coordinates: x, y, and z\n" + \
              "# vtu step: %d\n" % vtu_step + \
              "# time: %.4e\n" % _time
    for i in range(len(trench_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(trench_coords_x[i]) + " " + str(trench_coords_y[i]) + " " + str(trench_coords_z[i])
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
        filename = "trench_%05d.txt" % vtu_snapshot
    elif crust_only == 2:
        filename = "trench_l_%05d.txt" % vtu_snapshot
    else:
        filename = "trench_lu_%05d.txt" % vtu_snapshot
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions (upper boundary): %s" % fileout)

    # b. Bottom trench positions
    trench_b_coords_x, trench_b_coords_y, trench_b_coords_z = VtkP.ExportSlabInfoB()
    outputs = "# trench coordinates: x, y, and z\n" + \
              "# vtu step: %d\n" % vtu_step + \
              "# time: %.4e\n" % _time
    for i in range(len(trench_b_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(trench_b_coords_x[i]) + " " + str(trench_b_coords_y[i]) + " " + str(trench_b_coords_z[i])
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
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

    # Write the profile at the center of the slab
    center_profile_x, center_profile_y, center_profile_z = VtkP.ExportCenterProfile()
    outputs = "# center profile coordinates: x, y, and z\n" + \
              "# vtu step: %d\n" % vtu_step + \
              "# time: %.4e\n" % _time
    for i in range(len(center_profile_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(center_profile_x[i]) + " " + str(center_profile_y[i]) + " " + str(center_profile_z[i])
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
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

    # Write pinned depth information if required
    # a. Top pinned points
    pinD_coords_x, pinD_coords_y, pinD_coords_z = VtkP.ExportPinDInfo()
    outputs = "# pinned points at depth %.2f km coordinates: x, y, and z\n" % (pin_depth / 1e3) + \
              "# vtu step: %d\n" % vtu_step + \
              "# time: %.4e\n" % _time
    for i in range(len(pinD_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(pinD_coords_x[i]) + " " + str(pinD_coords_y[i]) + " " + str(pinD_coords_z[i])
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
        filename = "trench_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    elif crust_only == 2:
        filename = "trench_l_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    else:
        filename = "trench_lu_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions at depth %.2f km: %s" % (pin_depth / 1e3, fileout))

    # b. Bottom pinned points
    pinD_b_coords_x, pinD_b_coords_y, pinD_b_coords_z = VtkP.ExportPinDInfoB()
    outputs = "# pinned points at depth %.2f km coordinates: x, y, and z\n" % (pin_depth / 1e3) + \
              "# vtu step: %d\n" % vtu_step + \
              "# time: %.4e\n" % _time
    for i in range(len(pinD_b_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(pinD_b_coords_x[i]) + " " + str(pinD_b_coords_y[i]) + " " + str(pinD_b_coords_z[i])
    if crust_only == 1:
        # Mark the composition used for the envelope in the filename
        filename = "trench_b_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    elif crust_only == 2:
        filename = "trench_l_b_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    else:
        filename = "trench_lu_b_d%.2fkm_%05d.txt" % (pin_depth / 1e3, vtu_snapshot)
    fileout = os.path.join(output_path, filename)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions at depth %.2f km: %s" % (pin_depth / 1e3, fileout))
    end = time.time()
    print("Write trench positions, takes %.2f s" % (end - start))
    start = end

    # Write a slice at a fixed depth
    # todo_hv: Iterate over specified depths and export slab surface data
    for hv_depth in horizontal_velocity_depths:
        coords_x, coords_y, coords_z = VtkP.SlabSurfaceAtDepth(hv_depth)
        outputs = "# trench coordinates: x, y, and z\n" + \
                  "# vtu step: %d\n" % vtu_step + \
                  "# time: %.4e\n" % _time
        for i in range(len(coords_x)):
            if i > 0:
                outputs += "\n"
            outputs += str(coords_x[i]) + " " + str(coords_y[i]) + " " + str(coords_z[i])
        if crust_only == 1:
            # Mark the composition used for the envelope in the filename
            filename = "slab_surface_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth / 1e3)
        elif crust_only == 2:
            filename = "slab_surface_l_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth / 1e3)
        else:
            filename = "slab_surface_lu_%05d_d%.2fkm.txt" % (vtu_snapshot, hv_depth / 1e3)
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
    slab_envelop_interval_w = kwargs.get("slab_envelop_interval_w", 10e3)
    slab_envelop_interval_d = kwargs.get("slab_envelop_interval_d", 10e3)
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
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK(base_name, SlabMorphology, if_rewrite=if_rewrite, slab_envelop_interval_w=slab_envelop_interval_w,\
                                                slab_envelop_interval_d=slab_envelop_interval_d, slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only,\
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

def GetSlabDipAngle(case_dir, vtu_snapshot, Visit_Options, depth_lookup, depth_interval, **kwargs):
    '''
    Get slab dip angle at a given snapshot
    currently only works for slab center
    '''
    # options
    indent = kwargs.get("indent", 0)
    silence = kwargs.get("silence", False)
    crust_only = kwargs.get("crust_only", 0)
    
    # mark the composition used for the envelop in the filename
    if crust_only == 1:
        filename1 = "center_profile_%05d.txt" % vtu_snapshot
    elif crust_only == 2:
        filename1 = "center_profile_l_%05d.txt" % vtu_snapshot
    else:
        filename1 = "center_profile_lu_%05d.txt" % vtu_snapshot
    filein1 = os.path.join(case_dir, "vtk_outputs", filename1)
    if not os.path.isfile(filein1):
        if not silence:
            warnings.warn('%s%s: File %s is not found' % (" "*indent, Utilities.func_name(), filein1), UserWarning)
        return np.nan
    else:
        if not silence:
            print("%s%s: Found flle %s" % (" "*indent, Utilities.func_name(), filename1))
   
    # get the slab tip depth
    data1 = np.loadtxt(filein1)
    xs = data1[:, 0]
    zs = data1[:, 2]
    sfunc = interp1d(zs/1e6, xs/1e6, assume_sorted=True, fill_value="extrapolate")
  
    # interpret the x coordinate around the depth_lookup
    # and compute the dip angle
    # lhy11009: for some odd reason, the buildin interp function doesn't work
    Ro = Visit_Options.options['OUTER_RADIUS']
    geometry = Visit_Options.options['GEOMETRY']
    depth0 = depth_lookup - depth_interval
    z0 = Ro - depth0
    # x0 = np.interp(z0, zs, xs)
    # x0 = sfunc(z0/1e6) * 1e6
    x0 = 0.0
    for i in range(zs.size-1):
        if z0 < zs[i] and z0 > zs[i+1]:
            x0 = xs[i+1] * (z0 - zs[i]) / (zs[i+1] - zs[i]) + xs[i] * (z0 - zs[i+1]) / (zs[i] - zs[i+1])
    depth1 = depth_lookup
    z1 = Ro - depth1
    # x1 = np.interp(z1, zs, xs)
    # x1 = sfunc(z1/1e6) * 1e6
    x1 = 0.0
    for i in range(zs.size-1):
        if z1 < zs[i] and z1 > zs[i+1]:
            x1 = xs[i+1] * (z1 - zs[i]) / (zs[i+1] - zs[i]) + xs[i] * (z1 - zs[i+1]) / (zs[i] - zs[i+1])
    if not silence:
        print("%s(x0, z0): (%.4e, %.4e)" % (" "*indent, x0, z0))
        print("%s(x1, z1): (%.4e, %.4e)" % (" "*indent, x1, z1))
    dip = get_dip(x0, z0, x1, z1, geometry)

    # return 
    return dip


# todo_3d_chunk
def GetSlabDipAngle(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    silence = kwargs.get("silence", False)
    y_query = kwargs.get("y_query", 0.0)
    pin_depth = kwargs.get("pin_depth", 100e3)
    crust_only = kwargs.get("crust_only", 0)
    # is_chunk = geometry == "chunk"

    # initiation
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)
   
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

    return ts, depths, dips


def PlotSlabDipAngle(case_dir, time_interval_for_slab_morphology, **kwargs):
    '''
    plot trench position for a single stepb
    '''
    _color = kwargs.get("color", None)
    _label = kwargs.get("label", None)
    y_query = kwargs.get("y_query", 0.0)
    ax_twinx = kwargs.get("axis_twinx", None)
    # initiate plot
    ax = kwargs.get('axis', None)
    if ax == None:
        raise ValueError("Not implemented")
    ts, depths, dips = GetSlabDipAngle(case_dir, time_interval_for_slab_morphology, **kwargs)
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
    Visit_Options.Interpret(graphical_type="slice_center")
    trench_initial_position = Visit_Options.options['TRENCH_INITIAL']
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    _, available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology_outputs(time_interval=time_interval_for_slab_morphology)
    n_snapshots = len(available_pvtu_snapshots)
    if n_snapshots <= 1:
        raise FileNotFoundError("available snapshots are not found / only 1 snapshot for case %s" % case_dir)
    # derive the depth of the slab tip
    depths = []
    ts = []
    steps = []
    # n_snapshots - 1: some time the last step is not complete?
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
    try:
        t660 = sfunc(660e3)
        i660 = IndexByValue(ts, t660)
        step660 = steps[i660]
    except ValueError:
        t660 = np.nan
        i660 = np.nan
        step660 = np.nan
    try:
        t800 = sfunc(800e3)
        i800 = IndexByValue(ts, t800)
        step800 = steps[i800]
    except ValueError:
        t800 = np.nan
        i800 = np.nan
        step800 = np.nan
    try:
        t1000 = sfunc(1000e3)
        i1000 = IndexByValue(ts, t1000)
        step1000 = steps[i1000]
    except ValueError:
        t1000 = np.nan
        i1000 = np.nan
        step1000 = np.nan

    # construct the results
    results = {}
    results['t660'] = t660
    results['step660'] = step660

    results['t800'] = t800
    results['step800'] = step800
    
    results['t1000'] = t1000
    results['step1000'] = step1000

    return results


# todo_i3d
def Interpolate3dVtkCaseByParts(case_dir, vtu_snapshot, **kwargs):
    '''
    interpolate a 3d vtk output from a case, start from a single part
    '''
    assert(os.path.isdir(case_dir))
    # get options
    n0 = kwargs.get("n0", 800)
    n1 = kwargs.get("n1", 100)
    N = n0 * n1
    d_lateral = kwargs.get("d_lateral", 1.0)
    field = kwargs.get("field", "T")
    file_extension = kwargs.get("file_extension", 'txt')
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    resampled_data = np.zeros([target_points_np.shape[0], 1])
    resampled_data[:] = np.nan
    # rewrite the kwargs for sub-function 
    kwargs['file_extension'] = None # set to no output
    kwargs['target_points_np'] = target_points_np
    kwargs['indent'] = 4
    points = None; 
    point_data = None
    for part in range(16):
        print("Mask data from part %d" % part)
        points_part, point_data_part = Mask3dVtkCase(case_dir, vtu_snapshot, part, **kwargs)
        print("Get %d points from part %d" % (points_part.shape[0], part))
        if part == 0:
            points = points_part
            point_data = point_data_part
        else:
            points = np.concatenate([points, points_part], axis=0)
            point_data = np.concatenate((point_data, point_data_part), axis=0)
    print("Get %d points in total" % (points.shape[0]))
    print("Start interpolation")
    start = time.time()
    resampled_data = griddata(points, point_data, target_points_np, method='linear')
    end = time.time()
    print("End interpolation, take %.2f s" % (end-start))
#    for part in range(16):
#        print("Start part %d/16" % part)
#        kwargs['part'] = part
#        resampled_data_pt = Interpolate3dVtkCase(case_dir, vtu_snapshot, **kwargs)
#        mask = (~np.isnan(resampled_data_pt[:,0]))
#        resampled_data[mask] = resampled_data_pt[mask]
#        mask1 = (~np.isnan(resampled_data[:, 0]))
#        print("End part %d/16, Valid values in resampled data: %d/%d" % (part, mask1.sum(), N))
    # path to output
    odir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(odir):
        os.mkdir(odir)
    filebase = "center_slice"
    if d_lateral > 1.0e3:
        filebase = "lateral_slice_d%.2fkm" % (d_lateral)
    extension = "%05d.%s" % (vtu_snapshot, file_extension)
    fileout=os.path.join(odir, "%s-%s" % (filebase, extension))

    # initiate writer and output 
    if file_extension == 'txt':
        odata = np.concatenate((target_points_np, resampled_data), axis=1)
        np.savetxt(fileout, odata)
        # writer = vtk.vtkSimplePointsWriter()
        print("Saved file %s" % (fileout))
    elif file_extension == 'vtp':
        # output to vtp file
        o_poly_data = vtk.vtkPolyData()  # initiate poly data
        # new mesh
        target_cells_vtk = VtkPp.GetVtkCells2d(n0, n1)
        # insert points
        target_points_vtk = vtk.vtkPoints()
        for i in range(target_points_np.shape[0]):
          target_points_vtk.InsertNextPoint(target_points_np[i, 0], target_points_np[i, 1], target_points_np[i, 2])
        o_poly_data.SetPoints(target_points_vtk) # insert points
        o_poly_data.SetPolys(target_cells_vtk)
        # insert data
        resampled_field_vtk = numpy_to_vtk(resampled_data[:, 0], deep=1)
        resampled_field_vtk.SetName(field)
        o_poly_data.GetPointData().SetScalars(resampled_field_vtk)
        # write to file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileout)
        writer.SetInputData(o_poly_data)
        # writer.SetFileTypeToBinary()  # try this later to see if this works
        writer.Update()
        writer.Write()
        print("Saved file %s" % (fileout))


def Mask3dVtkCase(case_dir, vtu_snapshot, part, **kwargs):
    '''
    mask the 3d vtk data
    Inputs:
        kwargs:
            n0, n1 - number of points along the 1st and 3rd dimention
            interval - this determines the interval of the slices
            d_lateral - the lateral distance, along the 2nd dimention
                take a minimum value of 1.0 to assure successful slicing of the geometry
    '''
    # parse options
    interval = kwargs.get("interval", 10e3)
    field = kwargs.get("field", "T")
    file_extension = kwargs.get("file_extension", "vtp")
    n0 = kwargs.get("n0", 800)
    n1 = kwargs.get("n1", 100)
    indent = kwargs.get('indent', 0) # indentation when outputing to screen
    d_lateral = kwargs.get("d_lateral", 1.0)
    target_points_np = kwargs.get("target_points_np", None)
    print("%sMask data from part %d, get field %s" % (indent*" ", part, field))

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
    
    # parse file location, time and step
    # allow importing the solution for the step (.pvtu) or the solution for part of the stpe (.vtu)
    assert(type(part) == int)
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.%04d.vtu" % (vtu_snapshot, part))
    reader = vtk.vtkXMLUnstructuredGridReader()
    Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found" % filein)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)

    # parse from depth average files 
    Utilities.my_assert(_time != None, ValueError, "\"time\" is a requried input if \"ha_file\" is presented")

    # initiate vtk readers
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
    print("%sRead data takes %.2f s" % (" "*indent, end-start))
    print("%sBound of grid:" % (" "*indent), points.GetBounds())

    # get points and point data 
    points_np = vtk_to_numpy(points.GetData())
    field_np = vtk_to_numpy(point_data.GetArray(field))
    field_np = field_np.reshape((field_np.size, 1))

    if target_points_np is None:
        # make a target mesh is none is given
        target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    print("%sPoints in target: %d" % (" "*indent, target_points_np.shape[0]))
    # mask the data to the adjacency of the mesh
    mask=None
    if geometry == "box":
        mask = (points_np[:, 1] > (d_lateral - interval)) & (points_np[:, 1] < (d_lateral + interval))
    elif geometry == "chunk":
        mask = (points_np[:, 2] > (d_lateral - interval)) & (points_np[:, 2] < (d_lateral + interval))
    mask_points = points_np[mask, :]
    mask_data = field_np[mask, :]
    return mask_points, mask_data


def Interpolate3dVtkCase(case_dir, vtu_snapshot, **kwargs):
    '''
    interpolate a 3d vtk output from a case
    Inputs:
        kwargs:
            n0, n1 - number of points along the 1st and 3rd dimention
            interval - this determines the interval of the slices
            d_lateral - the lateral distance, along the 2nd dimention
                take a minimum value of 1.0 to assure successful slicing of the geometry
    '''
    interval = kwargs.get("interval", 10e3)
    file_extension = kwargs.get("file_extension", "vtp")
    n0 = kwargs.get("n0", 800)
    n1 = kwargs.get("n1", 100)
    indent = kwargs.get('indent', 0) # indentation when outputing to screen
    d_lateral = kwargs.get("d_lateral", 1.0)
    part = kwargs.get("part", None)
    target_points_np = kwargs.get("target_points_np", None)

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
    
    # load the class for processing the depth average output
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.Import(ha_file)

    # parse file location, time and step
    # allow importing the solution for the step (.pvtu) or the solution for part of the stpe (.vtu)
    filein = None
    reader = None
    if part is None:
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
        reader = vtk.vtkXMLPUnstructuredGridReader()
    else:
        assert(type(part) == int)
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.%04d.vtu" % (vtu_snapshot, part))
        reader = vtk.vtkXMLUnstructuredGridReader()
    Utilities.my_assert(os.path.isfile(filein), FileNotFoundError, "%s is not found" % filein)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)

    # parse from depth average files 
    Utilities.my_assert(_time != None, ValueError, "\"time\" is a requried input if \"ha_file\" is presented")
    Tref_func = DepthAverage.GetInterpolateFunc(_time, "temperature")
    density_ref_func = DepthAverage.GetInterpolateFunc(_time, "adiabatic_density")

    # initiate vtk readers
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
    print("%sRead data takes %.2f s" % (" "*indent, end-start))
    print("%sBound of grid:" % (" "*indent), points.GetBounds())

    # get points and point data 
    points_np = vtk_to_numpy(points.GetData())
    T_np = vtk_to_numpy(point_data.GetArray("T"))

    if target_points_np is None:
        # make a target mesh is none is given
        target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    print("%sPoints in target: %d" % (" "*indent, target_points_np.shape[0]))
    # mask the data to the adjacency of the mesh
    mask=None
    if geometry == "box":
        mask = (points_np[:, 1] > (d_lateral - interval)) & (points_np[:, 1] < (d_lateral + interval))
    elif geometry == "chunk":
        mask = (points_np[:, 2] > (d_lateral - interval)) & (points_np[:, 2] < (d_lateral + interval))

    # apply a mask before interpolate
    print("%sInterval = %.2e, Adjacent points number = %d" % (" "*indent, interval, points_np[mask].shape[0]))

    # Resample the data from the source mesh to the target mesh
    resampled_data = griddata(points_np[mask], T_np[mask], target_points_np, method='linear')
    if resampled_data.ndim == 1:
        # reshape to the right shape if there is only one field
        resampled_data = resampled_data.reshape(resampled_data.size, 1)
    mask = (~np.isnan(resampled_data[:, 0]))
    print("%sValid values in resampled data: %d" % (" "*indent, mask.sum()))

    # path to output
    odir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(odir):
        os.mkdir(odir)
    filebase = "center_slice"
    if d_lateral > 1.0e3:
        filebase = "lateral_slice_d%.2fkm" % (d_lateral)
    if part is None:
        extension = "%05d.%s" % (vtu_snapshot, file_extension)
    else:
        extension = "%05d.%04d.%s" % (vtu_snapshot, part, file_extension)
    fileout=os.path.join(odir, "%s-%s" % (filebase, extension))

    # initiate writer and output 
    if file_extension == 'txt':
        odata = np.concatenate((target_points_np, resampled_data), axis=1)
        np.savetxt(fileout, odata)
        # writer = vtk.vtkSimplePointsWriter()
        print("%sSaved file %s" % (" "*indent, fileout))
    elif file_extension == 'vtp':
        # output to vtp file
        o_poly_data = vtk.vtkPolyData()  # initiate poly data
        # new mesh
        target_cells_vtk = VtkPp.GetVtkCells2d(n0, n1)
        # insert points
        target_points_vtk = vtk.vtkPoints()
        for i in range(target_points_np.shape[0]):
          target_points_vtk.InsertNextPoint(target_points_np[i, 0], target_points_np[i, 1], target_points_np[i, 2])
        o_poly_data.SetPoints(target_points_vtk) # insert points
        o_poly_data.SetPolys(target_cells_vtk)
        # insert data
        resampled_T_vtk = numpy_to_vtk(resampled_data[:, 0], deep=1)
        resampled_T_vtk.SetName('T')
        o_poly_data.GetPointData().SetScalars(resampled_T_vtk)
        # write to file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileout)
        writer.SetInputData(o_poly_data)
        # writer.SetFileTypeToBinary()  # try this later to see if this works
        writer.Update()
        writer.Write()
        print("%sSaved file %s" % (" "*indent, fileout))
    else:
        pass
    
    return resampled_data


# todo_3d_chunk
def get_slab_dimensions_3(x, y, z, Ro, is_chunk):
    '''
    Derives the length along the three dimensions of a subducting slab.

    Inputs:
        x (float): x-coordinate of the slab point.
        y (float): y-coordinate of the slab point.
        z (float): z-coordinate of the slab point.
        Ro (float): Outer radius of the spherical domain.
        is_chunk (bool): Flag indicating whether the geometry is a spherical chunk.

    Returns:
        tuple: A tuple containing (r, w, l):
            - r (float): Radius or z-coordinate depending on whether the geometry is a chunk.
            - w (float): Width of the slab in the y-dimension, or converted width for chunk geometry.
            - l (float): Length of the slab in the x-dimension, or converted length for chunk geometry.
    
    Description:
        - For chunk geometries, converts Cartesian coordinates to spherical coordinates and calculates
          width and length using the outer radius Ro and spherical angles.
        - For non-chunk geometries, returns the z, x, and y coordinates directly as radius, length, and width.
    '''
    if is_chunk:
        # Convert Cartesian coordinates to spherical coordinates for chunk geometry
        r, th1, ph1 = Utilities.cart2sph(x, y, z)
        w = Ro * (np.pi / 2.0 - th1)  # Calculate width using the spherical angle th1
        l = Ro * ph1  # Calculate length using the spherical angle ph1
    else:
        # For non-chunk geometry, use Cartesian coordinates directly
        r = z
        l = x
        w = y

    return r, w, l


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
    parser.add_argument('-seix', '--slab_envelop_interval_w', type=float,
                        default=20e6,
                        help='Interval along y axis to sort out the trench locations')
    parser.add_argument('-seiz', '--slab_envelop_interval_d', type=float,
                        default=20e6,
                        help='Interval along z axis to sort out the trench locations')
    parser.add_argument('-ssc', '--slab_shallow_cutoff', type=float,
                        default=40e6,
                        help='Minimum depth along z axis to sort out the trench locations')
    parser.add_argument('-co', '--crust_only', type=int,
                        default=0,
                        help='If we only use the crustal composition to sort out the trench locations')
    parser.add_argument('-sptb', '--split_perturbation', type=int,
                        default=2,
                        help='his determines the number of slices the algorithm searches for a query point')
    parser.add_argument('-gs', '--geometry_spacing', nargs='+', type=int,
                        default=[10, 10, 10],
                        help='a list of values, spacing of the domain, this determines the number of slices the interpolation algorithm produce')
    parser.add_argument('-sr', '--slice_resolution', nargs='+', type=int,
                        default=[800, 300],
                        help='a list of values, mesh resolutions for interpolation')
    parser.add_argument('-sdl', '--slice_d_lateral', type=float,
                        default=1e3,
                        help='this is the distance of the slice measured on the 2nd dimension')
    parser.add_argument('-aacs', '--apply_additional_chunk_search', type=int,
                        default=1,
                        help='Apply additional search in interpolation if interpolatiion is in chunk geometry.')
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
            slab_envelop_interval_w=arg.slab_envelop_interval_w, slab_envelop_interval_d=arg.slab_envelop_interval_d,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only, include_force_balance=include_force_balance)
    elif _commend == 'cross_section_at_depth':
        SlabMorphology(arg.inputs, int(arg.vtu_snapshot), rewrite=1,\
            slab_envelop_interval_w=arg.slab_envelop_interval_w, slab_envelop_interval_d=arg.slab_envelop_interval_d,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only, horizontal_velocity_depths=[arg.depth])
    elif _commend == 'morph_case':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=0, time_interval=arg.time_interval,\
            slab_envelop_interval_w=arg.slab_envelop_interval_w, slab_envelop_interval_d=arg.slab_envelop_interval_d,\
            slab_shallow_cutoff=arg.slab_shallow_cutoff, crust_only=arg.crust_only)
    elif _commend == 'morph_case_parallel':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=1, time_interval=arg.time_interval,\
            slab_envelop_interval_w=arg.slab_envelop_interval_w, slab_envelop_interval_d=arg.slab_envelop_interval_d,\
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
    # todo_inter
    elif _commend == 'slice_3d_geometry_alpha':
        # Interpolate3dVtkCase(arg.inputs, arg.vtu_snapshot, interval=10e3, n0=800, n1=100, file_extension="vtp")
        # Interpolate3dVtkCase(arg.inputs, arg.vtu_snapshot, interval=1000e3, n0=800, n1=100, file_extension="txt", part=2)
        # Interpolate3dVtkCaseByParts(arg.inputs, arg.vtu_snapshot, interval=1000e3, n0=800, n1=100, file_extension="txt")
        Interpolate3dVtkCaseByParts(arg.inputs, arg.vtu_snapshot, interval=10e3, n0=800, n1=100, file_extension="vtp", field="viscosity")
        Interpolate3dVtkCaseByParts(arg.inputs, arg.vtu_snapshot, interval=10e3, n0=800, n1=100, d_lateral=750e3, file_extension="vtp", field="viscosity")
    elif _commend == 'slice_3d_geometry':
        fields = ["viscosity", "T"]
        mesh_options = {"type": "slice_2nd", "resolution": arg.slice_resolution, "d_lateral": arg.slice_d_lateral}
        print("mesh_options: ", mesh_options)
        print("spacing: ", arg.geometry_spacing)
        print("split_perturbation: ", arg.split_perturbation)
        print("apply_additional_chunk_search: ", arg.apply_additional_chunk_search)
        VtkPp.Interpolate3dVtkCaseBeta(arg.inputs, VISIT_OPTIONS, arg.vtu_snapshot, fields, mesh_options,\
             spacing=arg.geometry_spacing, split_perturbation=arg.split_perturbation, by_part=True,\
                 apply_additional_chunk_search=(arg.apply_additional_chunk_search==1))
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()