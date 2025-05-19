# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
"""
#### 3rd parties 
import numpy as np
import pandas as pd
import sys, os, argparse, json
# import json, re
# import pathlib
# import subprocess
import vtk
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from matplotlib import gridspec, cm
from matplotlib import colors as mcolors
from cmcrameri import cm as ccm
import shilofue.VtkPp as VtkPp
from shilofue.VtkPp import get_r
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from numpy import linalg as LA 
import multiprocessing
import warnings
import time
from scipy.interpolate import interp1d, griddata, UnivariateSpline
from scipy.integrate import quad
from scipy.optimize import minimize
#### self
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.PlotCombine import PLOT_COMBINE, PlotCombineExecute, PlotColorLabels, PC_OPT_BASE
# from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile, ParseFromDealiiInput, FindWBFeatures
from shilofue.Plot import LINEARPLOT
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

HaMaGeoLib_path = os.path.abspath("/home/lochy/ASPECT_PROJECT/HaMaGeoLib")
if HaMaGeoLib_path not in sys.path:
    sys.path.append(HaMaGeoLib_path)

from hamageolib.research.haoyuan_2d_subduction.legacy_tools import SlabTemperature, SlabTemperatureCase

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - Write slab forces output: \n\
\n\
        python -m shilofue.TwoDSubduction0.VtkPp analyze_slab\n\
            -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 -vss 105 -o foo.txt\n\
\n\
  - plot the slab forces: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_forces -i foo.txt -o foo.png\n\
\n\
  - Perform analysis on a single step for one case: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_case_step -i \n\
            /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0  -vs 100\n\
\n\
  - Slab morphology for a single step of a single case, by default, the slab envelops are outputed for debug usage and previous results are overwritten.\n\
        python -m shilofue.TwoDSubduction0.VtkPp morph_step -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 -vss 105 \n\
\n\
  - Slab morphology for one case for all the time steps\n\
        python -m shilofue.TwoDSubduction0.VtkPp morph_case -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT_cart2/eba_cdpt_cart_width80 -ti 0.1e6\n\
        note the -ti option would assign an interval (default is 0.5e6)\n\
\n\
  - Slab morphology for one case for all the time steps in parallel\n\
        python -m shilofue.TwoDSubduction0.VtkPp morph_case_parallel -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT_cart2/eba_cdpt_cart_width80 -ti 0.1e6\n\
        note the -ti option would assign an interval (default is 0.5e6)\n\
\n\
  - Plot the morphology of the slab: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_morph -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 \n\
\n\
  - Plot the morphology for publication: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_morph_publication -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT18/eba_cdpt_coh500_SA80.0_OA40.0_cd100.0_cd7.5_gr10\n\
\n\
  - Combine the morphology of a few cases. Note for this to work, a json file needs to be presented: \n\
        python -m shilofue.TwoDSubduction0.VtkPp combine_slab_morph -j /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT_peierls1/plot_combine/slab_morph.json\n\
\n\
  - Generature figure of shear zone geometry for a case at one step: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_shear_zone -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -vss 100\n\
\n\
  - Generature figure of shear zone geometry for a case, with a time interval and an end: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_shear_zone_case -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ti 2e6 -te 10e6\n\
\n\
        note the -ti option would assign an interval and the -te option will assign and end\n\
\n\
  - Generate a plot of crustal material over depth in each segment: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_material_time -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT4/eba_cdpt_SA80.0_OA40.0_CpEcl -t 1e6\n\
        the -t option would assign a model time to plot\n\
\n\
  - Generate a plot of temperature on the slab over depth in each segment: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_temperature -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -vss 100\n\
\n\
  - Generate the outputs of slab temperature\n\
        python -m shilofue.TwoDSubduction0.VtkPp slab_temperature_case -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140\n\
\n\
  - Plot the outputs of slab temperature\n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_temperature_case -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ts 10e6 -te 60e6\n\
\n\
  - Generate the outputs of mantle wedge temperature\n\
        python -m shilofue.TwoDSubduction0.VtkPp mantle_wedge_T_case -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ti 0.5e6\n\
\n\
  - Plot the mantle wedge temperature\n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_wedge_T -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ti 0.5e6\n\
        (note the interval of time assigned needs to be the same as the mantle_wedge_T_case command)\n\
\n\
  - Plot the thermal state of the slab\n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_thermal -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ti 0.5e6\n\
        (note the slab_morph.txt must be present)\n\
\n\
  - Plot the slab age intepreted from the thermal state of the trench\n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_age_from_T -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0_width140 -ti 0.5e6\n\
        (note the trench_T.txt files need to be present)\n\
        ")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
        slab_cells: cell id of internal points in the slab
        slab_envelop_cell_list0: cell id of slab envelop with smaller theta, slab bottom
        slab_envelop_cell_list1: cell id of slab envelop with bigger theta, slab surface
        trench: trench position, theta for a 'chunk' model and x for a 'box' model
        coord_100: position where slab is 100km deep, theta for a 'chunk' model and x for a 'box' model
        vsp : velocity of the subducting plate
        vov : velocity of the overiding plate
    '''
    def __init__(self, **kwargs):
        '''
        Initiation
        kwargs (dict):
            geometry - type of geometry
            Ro - outer radius
        '''
        VtkPp.VTKP.__init__(self, **kwargs)
        self.slab_cells = []
        self.crust_cells = []
        self.surface_cells = []
        self.slab_envelop_cell_list0 = []
        self.slab_envelop_cell_list1 = []
        self.cmb_envelop_cell_list = []
        self.sz_geometry = None
        self.slab_depth = None
        self.trench = None
        self.coord_100 = None 
        self.coord_200 = None # relate with the previous one
        self.dip_100 = None
        self.vsp = None
        self.vov = None
        self.coord_distant_200 = None
        self.v_distant_200 = None
        self.visc_distant_200 = None
        self.slab_shallow_cutoff = kwargs.get("slab_shallow_cutoff", 50e3)  # depth limit to slab
        self.slab_envelop_interval = kwargs.get("slab_envelop_interval", 5e3)
        self.velocitw_query_depth = 5e3  # depth to look up plate velocities
        self.velocitw_query_disl_to_trench = 500e3  # distance to trench to look up plate velocities
        self.shallow_trench = None
        default_gravity_file = os.path.join(Utilities.var_subs('${ASPECT_SOURCE_DIR}'),\
        "data", "gravity-model", "prem.txt") 
        gravity_file = kwargs.get('gravity_file', default_gravity_file)
        assert(os.path.isfile(gravity_file))
        self.ImportGravityData(gravity_file)

    # todo_shallow
    def PrepareSlabShallow(self, trench_initial, **kwargs):
        '''
        Prepares and identifies shallow and bottom points of a slab for trench analysis,
        exporting relevant data to a file if specified. Calculates original and corrected
        distances between shallow and bottom points based on threshold values.
    
        Parameters:
            trench_initial (float): Initial trench position in radians.
            **kwargs: Additional options, including trench_lookup_range and export_shallow_file.
                n_crust - number of crust in the model
    
        Returns:
            dict: Contains information on original and corrected points with distances.
        '''
        shallow_cutoff = 1e3 # m
        bottom_start = 7.5e3
        bottom_cutoff = 30e3
        pinned_field_value_threshold = 0.8
    
        trench_lookup_range = kwargs.get("trench_lookup_range", 10.0 * np.pi / 180.0)
        export_shallow_file = kwargs.get("export_shallow_file", None)
        n_crust = kwargs.get("n_crust", 1)
    
        print("%s started" % Utilities.func_name())
        start = time.time()
    
        # Sorts the shallow points and retrieves relevant data arrays for processing.
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        point_data = self.i_poly_data.GetPointData()

        pinned_field = vtk_to_numpy(point_data.GetArray("spcrust"))
        pinned_bottom_field = vtk_to_numpy(point_data.GetArray("spharz"))
    
        shallow_points_idx = []
        bottom_points_idx = []
        for idx in range(self.i_poly_data.GetNumberOfPoints()):
            x, y, _ = points[idx]
            r = None; th = None; ph = None
            if self.is_chunk:
                r, th, ph = Utilities.cart2sph(x,y,0.0)
            else:
                r = y
            d = self.Ro - r
            pinned_field_value = pinned_field[idx]
            pinned_bottom_field_value = pinned_bottom_field[idx]
    
            # Determines points based on depth and field thresholds, and location near trench.
            if d < shallow_cutoff and pinned_field_value > pinned_field_value_threshold and \
                ph > trench_initial - trench_lookup_range and ph < trench_initial + trench_lookup_range:
                shallow_points_idx.append(idx)
    
            if d > bottom_start and d < bottom_cutoff and pinned_bottom_field_value > pinned_field_value_threshold and \
                ph > trench_initial - trench_lookup_range and ph < trench_initial + trench_lookup_range:
                bottom_points_idx.append(idx)
    
        n_shallow_points = len(shallow_points_idx)
        n_bottom_points = len(bottom_points_idx)
    
        # Extracts coordinates of shallow and bottom points for further processing.
        shallow_points = np.zeros([n_shallow_points, 3])
        for i in range(n_shallow_points):
            idx = shallow_points_idx[i]
            shallow_points[i] = points[idx]
    
        bottom_points = np.zeros([n_bottom_points, 3])
        for i in range(n_bottom_points):
            idx = bottom_points_idx[i]
            bottom_points[i] = points[idx]
    
        end = time.time()
        print("\tSort shallow points, takes %.2f s" % (end - start))
        start = end
    
        # Exports shallow and bottom points to a specified file if export option is provided.
        if export_shallow_file is not None:
            export_shallow_outputs = np.zeros([n_shallow_points+n_bottom_points, 3])
            for i in range(n_shallow_points):
                idx = shallow_points_idx[i]
                export_shallow_outputs[i] = points[idx]
            for i in range(n_shallow_points, n_shallow_points+n_bottom_points):
                idx = bottom_points_idx[i-n_shallow_points]
                export_shallow_outputs[i] = points[idx]
            with open(export_shallow_file, 'w') as fout:
                np.savetxt(fout, export_shallow_outputs, header="X Y Z\n%d %d" % (n_shallow_points, n_bottom_points))
            print("\t%s: Write output file %s" % (Utilities.func_name(), export_shallow_file))
            end = time.time()
            print("\tWrite output file, takes %.2f s" % (end - start))
            start = end
    
        # Identifies furthest shallow point based on angle for distance calculations.
        Phi = np.zeros(n_shallow_points)
        for i in range(n_shallow_points):
            idx = shallow_points_idx[i]
            x, y, z = points[idx]
            r, th, ph = Utilities.cart2sph(x, y, z)
            Phi[i] = ph
    
        i_max = np.argmax(Phi)
        id_max = shallow_points_idx[i_max]
        phi_max = Phi[i_max]
    
        # Retrieves initial furthest point coordinates for calculating distances to bottom points.
        x_max, y_max, z_max = points[id_max]
        min_dist = float('inf')
        i_min, min_dist = minimum_distance_array(bottom_points, x_max, y_max, z_max)
        x_b_min, y_b_min, z_b_min = bottom_points[i_min]
    
        outputs = {}
        outputs["original"] = {"points": [x_max, y_max, z_max], "match points": [x_b_min, y_b_min, z_b_min], "distance": min_dist}
    
        end = time.time()
        print("\tCalculate original distance, takes %.2f s" % (end - start))
        start = end
    
        # Applies corrections to the furthest point based on incremental angle adjustments.
        if not self.is_chunk:
            # Ensures only chunk geometry is handled; raises exception otherwise.
            raise NotImplementedError()
    
        dphi = 0.001 * np.pi / 180.0
        dphi_increment = 0.001 * np.pi / 180.0
        phi_new = phi_max
        i_min = None
        while dphi < 1.0:
            phi = phi_max - dphi
            x_new, y_new, z_new = Utilities.ggr2cart(0.0, phi, self.Ro)
            i_min, min_dist = minimum_distance_array(bottom_points, x_new, y_new, z_new)
            if min_dist < 9e3:
                phi_new = phi
                break
            dphi += dphi_increment
 
        x_new, y_new, z_new = Utilities.ggr2cart(0.0, phi_new, self.Ro)
        x_b_min, y_b_min, z_b_min = bottom_points[i_min]
    
        outputs["corrected"] = {"points": [x_new, y_new, z_new], "match points": [x_b_min, y_b_min, z_b_min], "distance": min_dist}

        self.shallow_trench = [x_new, y_new, z_new]
    
        end = time.time()
        print("\tPerform correction, takes %.2f s" % (end - start))
        start = end
    
        return outputs

    def PrepareSlab(self, slab_field_names, **kwargs):
        '''
        Prepares slab composition by identifying slab and crustal regions in the model.
        
        Parameters:
            slab_field_names (list): Names of fields used to define slab composition.
            kwargs (dict): Optional keyword arguments:
                - prepare_cmb (str, optional): Field name for core-mantle boundary (CMB) preparation.
                - prepare_slab_distant_properties (optional): Additional properties for distant slab.
                - slab_threshold (float, optional): Threshold value for slab field to classify as slab material.
                - depth_lookup (float, optional): Depth for slab surface lookup, default is 100 km.
                - depth_distant_lookup (float, optional): Depth threshold for distant properties, default is 200 km.
        '''
        # Ensure cell centers are included in data
        assert(self.include_cell_center)

        # Initialize optional parameters from kwargs
        prepare_cmb = kwargs.get('prepare_cmb', None)
        prepare_slab_distant_properties = kwargs.get('prepare_slab_distant_properties', None)
        slab_threshold = kwargs.get('slab_threshold', 0.2)
        depth_lookup = kwargs.get("depth_lookup", 100e3)
        depth_distant_lookup = kwargs.get('depth_distant_lookup', 200e3)

        # Extract point and cell data arrays from VTK data structures
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        point_data = self.i_poly_data.GetPointData()
        cell_point_data = self.c_poly_data.GetPointData()
        
        # Get the primary slab field and initialize crust field if `prepare_cmb` is set
        slab_field = VtkPp.OperateDataArrays(cell_point_data, slab_field_names,\
        [0 for i in range(len(slab_field_names) - 1)])
        crust_field = None  # store the field of crust composition
        if prepare_cmb is not None:
            crust_field = vtk_to_numpy(cell_point_data.GetArray(prepare_cmb))
        
        # Identify cells based on composition and radius, storing minimum radius for slab depth
        min_r = self.Ro
        for i in range(self.i_poly_data.GetNumberOfCells()):
            cell = self.i_poly_data.GetCell(i)
            id_list = cell.GetPointIds()  # list of point ids in this cell
            x = centers[i][0]
            y = centers[i][1]
            r = get_r(x, y, self.geometry)
            slab = slab_field[i]
            if slab > slab_threshold and ((self.Ro - r) > self.slab_shallow_cutoff):
                self.slab_cells.append(i)
                if r < min_r:
                    min_r = r
        
        # If CMB preparation is requested, identify crustal cells similarly
        if prepare_cmb is not None:
            # cells of the crustal composition
            for i in range(self.i_poly_data.GetNumberOfCells()):
                cell = self.i_poly_data.GetCell(i)
                id_list = cell.GetPointIds()  # list of point ids in this cell
                x = centers[i][0]
                y = centers[i][1]
                r = get_r(x, y, self.geometry)
                crust = crust_field[i]
                if crust > slab_threshold and ((self.Ro - r) > self.slab_shallow_cutoff):
                    self.crust_cells.append(i)
        
        # Calculate slab depth based on minimum radius
        self.slab_depth = self.Ro - min_r

        # Group slab cells into envelop intervals based on radial depth
        total_en_interval = int((self.slab_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval + 1)
        slab_en_cell_lists = [ [] for i in range(total_en_interval) ]
        for id in self.slab_cells:
            x = centers[id][0]  # first, separate cells into intervals
            y = centers[id][1]
            r = get_r(x, y, self.geometry)
            id_en =  int(np.floor(
                                  (self.Ro - r - self.slab_shallow_cutoff)/
                                  self.slab_envelop_interval))# id in the envelop list
            slab_en_cell_lists[id_en].append(id)

        # Find angular boundaries (min and max theta) for each slab interval
        for id_en in range(len(slab_en_cell_lists)):
            theta_min = 0.0  # then, loop again by intervals to look for a
            theta_max = 0.0  # max theta and a min theta for each interval
            cell_list = slab_en_cell_lists[id_en]
            if len(cell_list) == 0:
                continue  # make sure we have some point
            is_first = True
            id_min = -1
            id_max = -1
            for id in cell_list:
                x = centers[id][0]
                y = centers[id][1]
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
            self.slab_envelop_cell_list0.append(id_min)  # first half of the envelop
            self.slab_envelop_cell_list1.append(id_max)  # second half of the envelop
        
        # Identify the trench position based on maximum angular position
        id_tr = self.slab_envelop_cell_list1[0] # point of the trench
        x_tr = centers[id_tr][0]  # first, separate cells into intervals
        y_tr = centers[id_tr][1]
        self.trench = get_theta(x_tr, y_tr, self.geometry)

        # Calculate dip angle at a specified depth lookup (e.g., 100 km)
        self.coord_100 = self.SlabSurfDepthLookup(depth_lookup)
        if self.geometry == "chunk":
            x100 = (self.Ro - depth_lookup) * np.cos(self.coord_100)
            y100 = (self.Ro - depth_lookup) * np.sin(self.coord_100)
        elif self.geometry == "box":
            x100 = self.coord_100
            y100 = self.Ro - depth_lookup
        r100 = get_r(x100, y100, self.geometry)
        theta100 = get_theta(x100, y100, self.geometry)
        self.dip_100 = get_dip(x_tr, y_tr, x100, y100, self.geometry)

        # If CMB is being prepared, repeat envelope grouping and interval checks for crust cells
        if prepare_cmb is not None:
            # get the crust envelop
            crust_en_cell_lists = [ [] for i in range(total_en_interval) ]
            for id in self.crust_cells:
                x = centers[id][0]  # first, separate cells into intervals
                y = centers[id][1]
                r = get_r(x, y, self.geometry)
                id_en =  int(np.floor(
                                    (self.Ro - r - self.slab_shallow_cutoff)/
                                    self.slab_envelop_interval))# id in the envelop list
                crust_en_cell_lists[id_en].append(id) 

            # Identify angular boundaries for crust cells to match slab envelope intervals
            for id_en in range(len(crust_en_cell_lists)):
                theta_min = 0.0  # then, loop again by intervals to look for a
                # theta_max = 0.0  # max theta and a min theta for each interval
                cell_list = crust_en_cell_lists[id_en]
                if len(cell_list) == 0:
                    if len(slab_en_cell_lists[id_en]) == 0:
                        pass
                    else:
                        # if there are points in the slab interval list
                        # I'll append some non-sense value here to make sure these 
                        # two have the same size
                        self.cmb_envelop_cell_list.append(-1)
                    continue  # make sure we have some point
                is_first = True
                id_min = -1
                # id_max = -1
                for id in cell_list:
                    x = centers[id][0]
                    y = centers[id][1]
                    theta = get_theta(x, y, self.geometry)  # cart
                    if is_first:
                        id_min = id
                        # id_max = id
                        theta_min = theta
                        # theta_max = theta
                        is_first = False
                    else:
                        if theta < theta_min:
                            id_min = id
                            theta_min = theta
                        # if theta > theta_max:
                        #     id_max = id
                        #    theta_max = theta
                self.cmb_envelop_cell_list.append(id_min)  # first half of the envelop

            # Ensure crust and slab envelopes have equal lengths 
            assert(len(self.cmb_envelop_cell_list)==len(self.slab_envelop_cell_list1))
    

    def PrepareSlabByDT(self, **kwargs):
        '''
        prepare slab composition by temperature difference to the reference adiabat
        Inputs:
            Tref_func: a function for the reference T profile.
        '''
        assert(self.include_cell_center)
        assert(self.Tref_func != None)
        slab_threshold = kwargs.get('slab_threshold', -100.0)
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        point_data = self.i_poly_data.GetPointData()
        cell_point_data = self.c_poly_data.GetPointData()
        # the temperature field
        T_field = vtk_to_numpy(cell_point_data.GetArray("T"))
        # add cells by composition
        min_r = self.Ro
        for i in range(self.i_poly_data.GetNumberOfCells()):
            cell = self.i_poly_data.GetCell(i)
            id_list = cell.GetPointIds()  # list of point ids in this cell
            x = centers[i][0]
            y = centers[i][1]
            r = get_r(x, y, self.geometry)
            Tref = self.Tref_func(self.Ro - r)
            T = T_field[i]
            if T - Tref < slab_threshold and ((self.Ro - r) > self.slab_shallow_cutoff):
                # note on the "<": slab internal is cold
                self.slab_cells.append(i)
                if r < min_r:
                    min_r = r
        self.slab_depth = self.Ro - min_r  # cart
        # get slab envelops
        total_en_interval = int((self.slab_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval + 1)
        slab_en_cell_lists = [ [] for i in range(total_en_interval) ]
        for id in self.slab_cells:
            x = centers[id][0]  # first, separate cells into intervals
            y = centers[id][1]
            r = get_r(x, y, self.geometry)
            id_en =  int(np.floor(
                                  (self.Ro - r - self.slab_shallow_cutoff)/
                                  self.slab_envelop_interval))# id in the envelop list
            slab_en_cell_lists[id_en].append(id)
        for id_en in range(len(slab_en_cell_lists)):
            theta_min = 0.0  # then, loop again by intervals to look for a
            theta_max = 0.0  # max theta and a min theta for each interval
            cell_list = slab_en_cell_lists[id_en]
            if len(cell_list) == 0:
                continue  # make sure we have some point
            is_first = True
            id_min = -1
            id_max = -1
            for id in cell_list:
                x = centers[id][0]
                y = centers[id][1]
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
            self.slab_envelop_cell_list0.append(id_min)  # first half of the envelop
            self.slab_envelop_cell_list1.append(id_max)  # second half of the envelop
        # trench
        id_tr = self.slab_envelop_cell_list1[0] # point of the trench
        x_tr = centers[id_tr][0]  # first, separate cells into intervals
        y_tr = centers[id_tr][1]
        self.trench = get_theta(x_tr, y_tr, self.geometry)
        # 100 km dip angle
        depth_lookup = 100e3
        self.coord_100 = self.SlabSurfDepthLookup(depth_lookup)
        if self.geometry == "chunk":
            x100 = (self.Ro - depth_lookup) * np.cos(self.coord_100)
            y100 = (self.Ro - depth_lookup) * np.sin(self.coord_100)
        elif self.geometry == "box":
            x100 = self.coord_100
            y100 = self.Ro - depth_lookup
        r100 = get_r(x100, y100, self.geometry)
        theta100 = get_theta(x100, y100, self.geometry)
        self.dip_100 = get_dip(x_tr, y_tr, x100, y100, self.geometry)
        pass

    def GetDipAtDepth(self, depth_lookup, depth_interval):
        # 100 km dip angle
        self.coord_0 = self.SlabSurfDepthLookup(depth_lookup-depth_interval)
        self.coord_1 = self.SlabSurfDepthLookup(depth_lookup)
        x_0, y_0, x_1, y_1 = None, None, None, None
        if self.geometry == "chunk":
            x_0 = (self.Ro - depth_lookup +depth_interval) * np.cos(self.coord_0)
            y_0 = (self.Ro - depth_lookup +depth_interval) * np.sin(self.coord_0)
            x_1 = (self.Ro - depth_lookup) * np.cos(self.coord_1)
            y_1 = (self.Ro - depth_lookup) * np.sin(self.coord_1)
        elif self.geometry == "box":
            x_0 = self.coord_0
            y_0 = self.Ro - depth_lookup + depth_interval
            x_1 = self.coord_1
            y_1 = self.Ro - depth_lookup
        r_0 = get_r(x_0, y_0, self.geometry)
        theta_0 = get_theta(x_0, y_0, self.geometry)
        r_1 = get_r(x_1, y_1, self.geometry)
        theta_1 = get_theta(x_1, y_1, self.geometry)
        dip = get_dip(x_0, y_0, x_1, y_1, self.geometry)
        print("x_0, y_0: ", x_0, y_0) # debug
        print("x_1, y_1: ", x_1, y_1)
        return dip

    def ExportOvAthenProfile(self, depth_distant_lookup, **kwargs):
        '''
        query a profile ending at depth_distant_lookup that is 5 deg to the
        Inputs:
            depth_distant_lookup - depth to place this query point
        '''
        project_velocity = kwargs.get("project_velocity", True)
        n_sample = kwargs.get("n_sample", 100)
        
        assert((self.coord_100 is not None) and (self.coord_100 is not None))

        # query the poly_data 
        query_grid = np.zeros((n_sample,2))
        v_distant_profile = np.zeros((n_sample, 2))
        depths = np.zeros(n_sample)
        value_fixings = np.zeros(n_sample)
        for i in range(n_sample):
            x_distant_200 = None; y_distant_200 = None
            depth = depth_distant_lookup * (1.0 * i / (n_sample-1.0))
            depths[i] = depth
            if self.geometry == "chunk":
                x_distant_200 = (self.Ro - depth) * np.cos(self.coord_100 + 5.0 * np.pi / 180.0)
                y_distant_200 = (self.Ro - depth) * np.sin(self.coord_100 + 5.0 * np.pi / 180.0)
            elif self.geometry == "box":
                x_distant_200 = self.coord_100 + 5.0 * np.pi / 180.0 * self.Ro
                y_distant_200 = self.Ro - depth
            query_grid[i, 0] = x_distant_200
            query_grid[i, 1] = y_distant_200
        # interpolate 
        query_poly_data = VtkPp.InterpolateGrid(self.i_poly_data, query_grid, quiet=True)
        query_vs = vtk_to_numpy(query_poly_data.GetPointData().GetArray('velocity'))
        query_viscs = vtk_to_numpy(query_poly_data.GetPointData().GetArray('viscosity'))
        # fix_missing_data
        # reason: interpolation fails upon chaning resolution
        # strategy: perturb by 0.1 degree at a time
        for i in range(n_sample):
            while (abs(query_viscs[i])) < 1e-6:
                # data missin
                value_fixings[i] += 1.0
                if self.geometry == "chunk":
                    # perturb by 0.1 degress
                    x_distant_200_1 = (self.Ro - depths[i]) * np.cos(self.coord_100 + (5.0 - 0.1*value_fixings[i]) * np.pi / 180.0)
                    y_distant_200_1 = (self.Ro - depths[i]) * np.sin(self.coord_100 + (5.0 - 0.1*value_fixings[i]) * np.pi / 180.0)
                    x_distant_200_2 = (self.Ro - depths[i]) * np.cos(self.coord_100 + (5.0 + 0.1*value_fixings[i]) * np.pi / 180.0)
                    y_distant_200_2 = (self.Ro - depths[i]) * np.sin(self.coord_100 + (5.0 + 0.1*value_fixings[i]) * np.pi / 180.0)
                elif self.geometry == "box":
                    x_distant_200_1 = self.coord_100 + (5.0 - 0.1*value_fixings[i]) * np.pi / 180.0 * self.Ro
                    y_distant_200_1 = self.Ro - depths[i]
                    x_distant_200_2 = self.coord_100 + (5.0 + 0.1*value_fixings[i]) * np.pi / 180.0 * self.Ro
                    y_distant_200_2 = self.Ro - depths[i]
                query_grid1 = np.zeros((1, 2))
                query_grid1[0, 0] = x_distant_200_1
                query_grid1[0, 1] = y_distant_200_1
                query_grid2 = np.zeros((1, 2))
                query_grid1[0, 0] = x_distant_200_2
                query_grid1[0, 1] = y_distant_200_2
                query_poly_data1 = VtkPp.InterpolateGrid(self.i_poly_data, query_grid1, quiet=True)
                query_poly_data2 = VtkPp.InterpolateGrid(self.i_poly_data, query_grid2, quiet=True)
                query_v1 = vtk_to_numpy(query_poly_data1.GetPointData().GetArray('velocity'))
                query_visc1 = vtk_to_numpy(query_poly_data1.GetPointData().GetArray('viscosity'))
                query_v2 = vtk_to_numpy(query_poly_data2.GetPointData().GetArray('velocity'))
                query_visc2 = vtk_to_numpy(query_poly_data2.GetPointData().GetArray('viscosity'))
                if (abs(query_visc1)) > 1e-6:
                    # fixed
                    query_vs[i, :] = query_v1
                    query_viscs[i] = query_visc1
                    query_grid[i, :] = query_grid1
                    break
                elif (abs(query_visc2)) > 1e-6:
                    query_vs[i, :] = query_v2
                    query_viscs[i] = query_visc2
                    query_grid[i, :] = query_grid2
                    break
        # project the velocity if needed and get the viscosity
        if project_velocity:
            # project velocity to theta direction in a spherical geometry
            v_distant_profile[:, 0], v_distant_profile[:, 1] = VtkPp.ProjectVelocity(x_distant_200, y_distant_200, query_vs, self.geometry)
        else:
            v_distant_profile[:, 0], v_distant_profile[:, 1] = query_vs[:, 0], query_vs[:, 1]
        
        return query_grid, depths, v_distant_profile, query_viscs, value_fixings
    
    
    def ExportSlabInfo(self):
        '''
        Output slab information
        '''
        return self.trench, self.slab_depth, self.dip_100

    def ExportVelocity(self, **kwargs):
        '''
        Output sp and ov plate velocity
        Inputs:
            kwargs:
                project_velocity - whether the velocity is projected to the tangential direction
        '''
        project_velocity = kwargs.get('project_velocity', False)
        assert(self.trench is not None)
        if self.geometry == "chunk":
            r_sp_query = self.Ro - self.velocitw_query_depth
            # theta_sp_query = self.trench - self.velocitw_query_disl_to_trench / self.Ro
            theta_sp_query = self.trench / 2.0
            r_ov_query = self.Ro - self.velocitw_query_depth
            # theta_ov_query = self.trench + self.velocitw_query_disl_to_trench / self.Ro
            theta_ov_query = (self.trench + self.Xmax) / 2.0
            x_sp_query = r_sp_query * np.cos(theta_sp_query)
            y_sp_query = r_sp_query * np.sin(theta_sp_query)
            x_ov_query = r_ov_query * np.cos(theta_ov_query)
            y_ov_query = r_ov_query * np.sin(theta_ov_query)
        elif self.geometry == "box":
            # x_sp_query = self.trench - self.velocitw_query_disl_to_trench
            x_sp_query = self.trench / 2.0
            y_sp_query = self.Ro - self.velocitw_query_depth
            # x_ov_query = self.trench + self.velocitw_query_disl_to_trench
            x_ov_query = (self.trench + self.Xmax) / 2.0
            y_ov_query = self.Ro - self.velocitw_query_depth
        query_grid = np.zeros((2,2))
        query_grid[0, 0] = x_sp_query
        query_grid[0, 1] = y_sp_query
        query_grid[1, 0] = x_ov_query
        query_grid[1, 1] = y_ov_query
        query_poly_data = VtkPp.InterpolateGrid(self.i_poly_data, query_grid, quiet=True)
        query_vs = vtk_to_numpy(query_poly_data.GetPointData().GetArray('velocity'))
        if project_velocity:
            # project velocity to theta direction in a spherical geometry
            self.vsp, _ = VtkPp.ProjectVelocity(x_sp_query, y_sp_query, query_vs[0, :], self.geometry)
            self.vov, _ = VtkPp.ProjectVelocity(x_ov_query, y_ov_query, query_vs[1, :], self.geometry)
        else:
            self.vsp = query_vs[0, :]
            self.vov = query_vs[1, :]
        return self.vsp, self.vov

    def SlabSurfDepthLookup(self, depth_lkp):
        '''
        Get point from the surface of the slab by depth
        '''
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        assert(len(self.slab_envelop_cell_list1) > 0)
        assert(depth_lkp < self.slab_depth)
        is_first = True
        coord_last = 0.0
        depth_last = 0.0
        for id in self.slab_envelop_cell_list1:
            x = centers[id][0]
            y = centers[id][1]
            r = get_r(x, y, self.geometry)
            coord = get_theta(x, y, self.geometry)
            depth = self.Ro - r
            if depth_last < depth_lkp and depth_lkp <= depth:
                coord_lkp = coord * (depth_lkp - depth_last) / (depth - depth_last) +\
                            coord_last * (depth_lkp - depth) / (depth_last - depth)
                break
            coord_last = coord
            depth_last = depth
        return coord_lkp

    
    def ExportSlabInternal(self, output_xy=False):
        '''
        export slab internal points
        '''
        cell_source = vtk.vtkExtractCells()
        cell_source.SetInputData(self.i_poly_data)
        cell_source.SetCellList(VtkPp.NpIntToIdList(self.slab_cells))
        cell_source.Update()
        slab_cell_grid = cell_source.GetOutput()
        if output_xy:
            coords = vtk_to_numpy(slab_cell_grid.GetPoints().GetData())
            return coords
        else:
            return slab_cell_grid
    
    def ExportSlabEnvelopCoord(self, **kwargs):
        '''
        export slab envelop envelops,
        outputs:
            coordinates in slab envelop
        '''
        assert (len(self.slab_envelop_cell_list0) > 0 and\
            len(self.slab_envelop_cell_list1) > 0)  # assert we have slab internels
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        slab_envelop0 = []
        slab_envelop1 = []
        # envelop 0
        xs = []
        ys = []
        for id in self.slab_envelop_cell_list0:
            x = centers[id][0]
            y = centers[id][1]
            xs.append(x)
            ys.append(y)
        slab_envelop0 = np.array([xs, ys])
        # envelop 1
        xs = []
        ys = []
        for id in self.slab_envelop_cell_list1:
            x = centers[id][0]
            y = centers[id][1]
            xs.append(x)
            ys.append(y)
        slab_envelop1 = np.array([xs, ys])
        return slab_envelop0.T, slab_envelop1.T

    def ExportSlabCmbCoord(self, **kwargs):
        '''
        export slab core-mantle boundary envelop
        returns:
            coordinates of points on the core-mantle boundary
        '''
        assert (len(self.cmb_envelop_cell_list) > 0)
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        xs = []
        ys = []
        for id in self.cmb_envelop_cell_list:
            x = -1e31  # these initial value are none sense
            y = -1e31  # but if there are negative ids, just maintain these values
            if id > 0:
                x = centers[id][0]
                y = centers[id][1]
            xs.append(x)
            ys.append(y)
        cmb_envelop = np.array([xs, ys])
        return cmb_envelop.T

    def ExportWedgeT(self, **kwargs):
        '''
        export the temperature in the mantle wedge
        Inputs:
            kwargs:
                fileout - path for output temperature, if this is None,
                        then no output is generated
        '''
        fileout = kwargs.get('fileout', None)
        depth_lookup = 100e3
        min_depth = 0.0
        max_depth = 100e3
        n_points = 100  # points for vtk interpolation
        o_n_points = 200  # points for output
        depths = np.linspace(0.0, max_depth, n_points)
        o_depths = np.linspace(0.0, max_depth, o_n_points)
        o_Ts = np.zeros(o_n_points)
        # look up for the point on the slab that is 100 km deep
        coord_lookup = self.SlabSurfDepthLookup(depth_lookup)
        vProfile2D = self.VerticalProfile2D((self.Ro - max_depth, self.Ro),\
                                            coord_lookup, n_points, fix_point_value=True)
        Tfunc = vProfile2D.GetFunction('T')
        outputs = "# 1: x (m)\n# 2: y (m)\n# 3: T\n"
        is_first = True
        for i in range(o_n_points):
            depth = o_depths[i]
            if self.geometry == "chunk":
                radius = self.Ro - depth
                x = radius * np.cos(coord_lookup)
                y = radius * np.sin(coord_lookup)
            elif self.geometry == "box":
                x = coord_lookup
                y = self.Ro - depth
            T = Tfunc(self.Ro - depth)
            o_Ts[i] = T
            if is_first:
                is_first = False
            else:
                outputs += "\n"
            outputs += "%.4e %.4e %.4e" % (x, y, T)
        if fileout is not None:
            with open(fileout, 'w') as fout:
                fout.write(outputs)
            print("%s: write file %s" % (Utilities.func_name(), fileout))  # screen output
        return o_depths, o_Ts

    def ExportTrenchT(self, **kwargs):
        '''
        export the temperature in the mantle wedge
        Inputs:
            trench: coordinate of trench position (theta or x)
            kwargs:
                fileout - path for output temperature, if this is None,
                        then no output is generated
        '''
        assert(self.trench is not None)  # assert that the trench posiiton is processed
        fileout = kwargs.get('fileout', None)
        distance_to_trench = 200e3
        # the point to look up needs to be on the the subducting plate
        if self.geometry == 'box':
            lookup = self.trench - distance_to_trench
        elif self.geometry == 'chunk':
            lookup = self.trench - distance_to_trench / self.Ro
        depth_lookup = 100e3
        min_depth = 0.0
        max_depth = 100e3
        n_points = 100  # points for vtk interpolation
        o_n_points = 200  # points for output
        depths = np.linspace(0.0, max_depth, n_points)
        o_depths = np.linspace(0.0, max_depth, o_n_points)
        o_Ts = np.zeros(o_n_points)
        # look up for the point on the slab that is 100 km deep
        vProfile2D = self.VerticalProfile2D((self.Ro - max_depth, self.Ro),\
                                            lookup, n_points, fix_point_value=True)
        Tfunc = vProfile2D.GetFunction('T')
        outputs = "# 1: x (m)\n# 2: y (m)\n# 3: depth (m)\n# 4: T (K)\n"
        is_first = True
        for i in range(o_n_points):
            depth = o_depths[i]
            if self.geometry == "chunk":
                radius = self.Ro - depth
                x = radius * np.cos(lookup)
                y = radius * np.sin(lookup)
            elif self.geometry == "box":
                x = lookup
                y = self.Ro - depth
            T = Tfunc(self.Ro - depth)
            o_Ts[i] = T
            if is_first:
                is_first = False
            else:
                outputs += "\n"
            outputs += "%.4e %.4e %.4e %.4e" % (x, y, depth, T)
        if fileout is not None:
            with open(fileout, 'w') as fout:
                fout.write(outputs)
            print("%s: write file %s" % (Utilities.func_name(), fileout))  # screen output
        return o_depths, o_Ts

    def SlabBuoyancy(self, v_profile, depth_increment):
        '''
        Compute the slab buoyancy
        Inputs:
            v_profile: vertical profile containing the reference profile
        Outputs:
            total_buoyancy: the total buoyancy forces in N/m
            b_profile: the depths and buoyancies in np array, 
                    the depths serve as the center in [depth - depth_increment/2.0, depth + depth_increment/2.0]
                    and the buoyancies contain a correspondent value for each range

        '''
        grav_acc = 10.0
        assert(self.include_cell_center)
        assert(len(self.slab_cells) > 0)
        n_depth = int(np.ceil(self.slab_depth / depth_increment))
        buoyancies = np.zeros(n_depth) # save the values of buoyancy with ranges of depths
        depths = []  # construct depths, each serve as the center in [depth - depth_increment/2.0, depth + depth_increment/2.0]
        for i in range(n_depth):
            depth = (i + 0.5) * depth_increment
            depths.append(depth)
        depths = np.array(depths)
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())  # note these are data mapped to cell center
        density_data = vtk_to_numpy(self.c_poly_data.GetPointData().GetArray('density'))
        density_ref_func = v_profile.GetFunction('density')
        # now compute the slab buoyancy
        total_buoyancy = 0.0
        for i in self.slab_cells:
            x = centers[i][0]
            y = centers[i][1]
            r = get_r(x, y, self.geometry)
            i_r = int(np.floor((self.Ro - r) / depth_increment))
            density = density_data[i]
            density_ref = density_ref_func(r)
            cell_size = self.cell_sizes[i]  # temp
            buoyancy = - grav_acc * (density - density_ref) * cell_size  # gravity
            buoyancies[i_r] += buoyancy
            total_buoyancy += buoyancy
        b_profile = np.zeros((n_depth, 2))
        b_profile[:, 0] = depths
        b_profile[:, 1] = buoyancies
        return total_buoyancy, b_profile

    def FindMDD(self, **kwargs):
        '''
        find the mechanical decoupling depth from the velocity field
        '''
        dx0 = kwargs.get('dx0', 10e3)
        dx1 = kwargs.get('dx1', 10e3)
        tolerance = kwargs.get('tolerance', 0.05)
        indent = kwargs.get("indent", 0)  # indentation for outputs
        debug = kwargs.get("debug", False) # output debug messages
        slab_envelop0, slab_envelop1 = self.ExportSlabEnvelopCoord()
        query_grid = np.zeros((2,2))
        mdd = -1.0
        for i in range(slab_envelop1.shape[0]):
            x = slab_envelop1[i, 0]
            y = slab_envelop1[i, 1]
            r = get_r(x, y, self.geometry) 
            theta = get_theta(x, y, self.geometry)
            if self.geometry == "chunk":
                x0 = r * np.cos(theta - dx0/self.Ro) 
                y0 = r * np.sin(theta - dx0/self.Ro)
                x1 = r * np.cos(theta + dx1/self.Ro)
                y1 = r * np.sin(theta + dx1/self.Ro)
            elif self.geometry == "box":
                x0 = x - dx0
                y0 = y
                x1 = x + dx1
                y1 = y
            query_grid[0, 0] = x0
            query_grid[0, 1] = y0
            query_grid[1, 0] = x1
            query_grid[1, 1] = y1
            query_poly_data = VtkPp.InterpolateGrid(self.i_poly_data, query_grid, quiet=True)
            query_vs = vtk_to_numpy(query_poly_data.GetPointData().GetArray('velocity'))
            vi = query_vs[0, :]
            vi_mag = (query_vs[0, 0]**2.0 + query_vs[0, 1]**2.0)**0.5
            vi_theta = np.arctan2(query_vs[0,0], query_vs[0, 1])
            vo = query_vs[1, :]
            vo_mag = (query_vs[1, 0]**2.0 + query_vs[1, 1]**2.0)**0.5
            vo_theta = np.arctan2(query_vs[1,0], query_vs[1, 1])
            depth = (self.Ro - (x**2.0 + y**2.0)**0.5)
            if debug:
                print("%sx: %.4e, y: %.4e, depth: %.4e, vi: [%.4e, %.4e] (mag = %.4e, theta = %.4e), vo: [%.4e, %.4e] (mag = %.4e, theta = %.4e)"\
                % (indent*" ", x, y, depth, vi[0], vi[1], vi_mag, vi_theta, vo[0], vo[1], vo_mag, vo_theta))
            if (abs((vo_mag - vi_mag)/vi_mag) < tolerance) and (abs((vo_theta - vi_theta)/vi_theta) < tolerance):
                mdd = depth
                break
        print("%sfindmdd_tolerance = %.4e, dx0 = %.4e, dx1 = %.4e" % (indent*" ", tolerance, dx0, dx1))
        if mdd > 0.0:
            print("%smdd = %.4e m" % (indent*" ", mdd))
            return mdd
        else:
            raise ValueError("FindMDD: a valid MDD has not been found, please considering changing the tolerance")
        pass

    def PrepareSZ(self, fileout, **kwargs):
        '''
        Get the geometry of the shear zone
        '''
        contour_poly_data = VtkPp.ExportContour(self.i_poly_data, 'spcrust', 0.9, fileout=fileout)
        contour_data = vtk_to_numpy(contour_poly_data.GetPoints().GetData())
        Dsz = kwargs.get("Dsz", 7.5e3)
        max_extent = 50e3
        min_extent = 2e3
        min_depth = 10e3
        max_depth = 200e3
        depth_interval = 2e3
        inf = 1e31 # define a vary big value
        small_depth_variation = 500 # a negalectable variation in depth
        theta_subgroup = np.pi / 10.0 # 1 degree
        n_group = int(np.floor(max_depth / depth_interval) + 1.0)
        idx_groups = [[] for i in range(n_group)]
        sz_depths  = []  # initial arrays
        sz_theta_mins = []
        sz_theta_maxs = []
        # step 1: a. assign points in groups, with an interval in depth
        #         b. append points on the surface
        for i in range(contour_data.shape[0]):
            x = contour_data[i, 0]
            y = contour_data[i, 1]
            r = get_r(x, y, self.geometry)
            # theta = get_theta(x, y, self.geometry) 
            depth = self.Ro - r
            if depth < max_depth and depth > min_depth:
                # assign points in groups, with an interval in depth
                i_group = int(np.floor(depth / depth_interval))
                idx_groups[i_group].append(i)
        # step 2: look for the points located near to a query point in depth
        for i_group in range(n_group-1):
            if len(idx_groups[i_group])==0:
                # skip groups that has no points
                continue
            query_depth = depth_interval * (i_group + 0.5)
            theta_min = np.pi
            theta_max = 0.0
            for i in (idx_groups[i_group]):
                # find the minimum and maximum theta in this group
                x = contour_data[i, 0]
                y = contour_data[i, 1]
                theta = get_theta(x, y, self.geometry) 
                if theta < theta_min and True:
                    theta_min = theta
                if theta > theta_max and True:
                    theta_max = theta
            extent = 0.0
            if self.geometry == "box":
                extent = theta_max - theta_min
            elif self.geometry == "chunk":
                extent = (self.Ro - depth) * (theta_max - theta_min)
            else:
                raise ValueError("Geometry must be either box or chunk")
            if extent < max_extent and extent > min_extent:
                # select the reasonable values of layer thickness
                sz_depths.append(query_depth)
                sz_theta_mins.append(theta_min)
                sz_theta_maxs.append(theta_max)
        # loop again, append some points
        sz_theta_min_start = sz_theta_mins[0]
        sz_depths_app = []
        sz_theta_mins_app = [] 
        sz_theta_maxs_app = []
        theta_max_moho = 0.0
        depth_at_theta_max_moho = 0.0
        found = False
        for i in range(contour_data.shape[0]):
            # find the points on the subducting plate moho, and this point has to be below the subducting
            # plate and near the trench
            x = contour_data[i, 0]
            y = contour_data[i, 1]
            r = get_r(x, y, self.geometry)
            theta = get_theta(x, y, self.geometry) 
            depth = self.Ro - r
            if abs(depth - Dsz) < small_depth_variation and theta < sz_theta_min_start and (theta-sz_theta_min_start)/sz_theta_min_start > -0.1:
                # print("theta: ",theta)
                found = True
                if (theta > theta_max_moho):
                    theta_max_moho = theta
                    depth_at_theta_max_moho = depth
        Utilities.my_assert(found, ValueError, "PrepareSZ: No point found on moho near the trench.")
        sz_depths_surf_shallow = []
        sz_theta_mins_surf_shallow = []
        sz_theta_maxs_surf_shallow = []
        for i in range(contour_data.shape[0]):
            # find the points on the subducting surface and shallower than the initial depth
            x = contour_data[i, 0]
            y = contour_data[i, 1]
            r = get_r(x, y, self.geometry)
            theta = get_theta(x, y, self.geometry) 
            depth = self.Ro - r
            if depth > 0.0 and depth < depth_at_theta_max_moho and theta > theta_max_moho:
                pass
                sz_depths_surf_shallow.append(depth)
                sz_theta_mins_surf_shallow.append(-inf)
                sz_theta_maxs_surf_shallow.append(theta)
        sz_depths_app += sz_depths_surf_shallow # append the points on the subducting surface and shallower than the initial depth
        sz_theta_mins_app += sz_theta_mins_surf_shallow
        sz_theta_maxs_app += sz_theta_maxs_surf_shallow
        sz_depths_app.append(0.0) # append two points, one on the surface, another on the moho
        sz_theta_mins_app.append(-inf)
        sz_theta_maxs_app.append(theta_max_moho)
        sz_depths_app.append(depth_at_theta_max_moho)
        sz_theta_mins_app.append(theta_max_moho)
        sz_theta_maxs_app.append(inf)
        sz_depths = sz_depths_app + sz_depths # append to the front
        sz_theta_mins = sz_theta_mins_app + sz_theta_mins
        sz_theta_maxs = sz_theta_maxs_app + sz_theta_maxs
        self.sz_geometry = np.array([sz_depths, sz_theta_mins, sz_theta_maxs]).transpose()
        # write output file
        header = "# 1: depth (m)\n# 2: theta_min\n# 3: theta_max\n"
        with open(fileout, 'w') as fout:
            fout.write(header)
            np.savetxt(fout, self.sz_geometry)
        print("%s: file output: %s" % (Utilities.func_name(), fileout))
                
####
# Utilities functions
####


def get_theta(x, y, geometry):
    '''
    Get theta (the second coordinate)
    Inputs:
        x - x coordinate
        y - y coordinate
        geometry - 'chunk' or 'box'
    '''
    if geometry == 'chunk':
        theta = np.arctan2(y, x)  # cart
    elif geometry == 'box':
        theta = x
    else:
        raise ValueError("not implemented")
    return theta


def get_dip(x0, y0, x1, y1, geometry):
    '''
    Get dip angle
    Inputs:
        x0, y0: coordinates of the first point
        x1, y1: coordinates of the second point
        geometry - 'chunk' or 'box'
    '''
    if geometry == 'chunk':
        # here, for the 2nd dimension, we need something multiple the change in theta,
        # and I pick (r1 + r0)/2.0 for it.
        theta0 = np.arctan2(y0, x0)  # cart
        theta1 = np.arctan2(y1, x1)  # cart
        dtheta = theta1 - theta0
        r0 = (x0*x0 + y0*y0)**0.5
        r1 = (x1*x1 + y1*y1)**0.5
        # dip = np.arctan2(r0-r1*np.cos(dtheta), r1*np.sin(dtheta))
        dip = np.arctan2(r0-r1, (r1 + r0)/2.0*dtheta)
    elif geometry == 'box':
        dip = np.arctan2(-(y1-y0), (x1-x0))
    else:
        raise ValueError("not implemented")
    return dip


class T_FIT_FUNC():
    '''Residual functions'''
    def __init__(self, depth_to_fit, T_to_fit, **kwargs):
        '''
        Inputs:
            depth_to_fit - an array of depth
            T_to_fit - an array of temperature to fit
            kwargs:
                potential_temperature - mantle potential temperature

        '''
        self.depth_to_fit = depth_to_fit
        self.T_fit = T_to_fit
        self.potential_temperature = kwargs.get('potential_temperature', 1673.0)

    def PlateModel(self, xs):
        '''
        Inputs:
            xs - an array of non-dimensioned variables
        '''
        Ts = []
        seconds_in_yr = 3600.0 * 24 * 365
        age = xs[0] * 40 * seconds_in_yr * 1e6
        for depth in self.depth_to_fit:
            Ts.append(plate_model_temperature(depth, age=age, potential_temperature=self.potential_temperature))
        return np.linalg.norm(Ts - self.T_fit, 2.0)


def plate_model_temperature(depth, **kwargs):
    '''
    Use plate model to compute the temperature, migrated from the world builder
    in order to generate consistent result with the world builder.
    Inputs:
        x - distance to ridge, in meter
        depth - depth under the surface
        kwargs:
            plate_velocity - velocity of the plate, if this value is given, then the
                            distance to the ridge is also needed.
            distance_ridge - distance to the ridge
            age - age of the plate
    '''
    plate_velocity = kwargs.get('plate_velocity', None)
    if plate_velocity is None:
        age = kwargs['age']
    else:
        distance_ridge = kwargs['distance_ridge']
        age = distance_ridge / plate_velocity
    max_depth = kwargs.get('max_depth', 150e3)
    potential_mantle_temperature = kwargs.get('potential_temperature', 1673.0)
    top_temperature = kwargs.get('top_temperature', 273.0)
    thermal_diffusivity = 1e-6
    thermal_expansion_coefficient = 3e-5
    gravity_norm = 10.0
    specific_heat = 1250.0
    sommation_number = 100 # same as in the World Builder
    bottom_temperature_local = potential_mantle_temperature *\
                            np.exp(((thermal_expansion_coefficient* gravity_norm) / specific_heat) * depth)

    temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth)

    for i in range(1, sommation_number+1):
        # suming over the "sommation_number"
        # use a spreading ridge around the left corner and a constant age around the right corner 
        if plate_velocity is not None:
            temperature = temperature + (bottom_temperature_local - top_temperature) *\
                        ((2.0 / (i * np.pi)) * np.sin((i * np.pi * depth) / max_depth) *\
                         np.exp((((plate_velocity * max_depth)/(2 * thermal_diffusivity)) -\
                                   np.sqrt(((plate_velocity*plate_velocity*max_depth*max_depth) /\
                                              (4*thermal_diffusivity*thermal_diffusivity)) + i * i * np.pi * np.pi)) *\
                                  ((plate_velocity * age) / max_depth)))
        else:
            temperature = temperature + (bottom_temperature_local - top_temperature) *\
                        ((2.0 / (i * np.pi)) * np.sin((i * np.pi * depth) / max_depth) *\
                         np.exp(-1.0 * i * i * np.pi * np.pi * thermal_diffusivity * age / (max_depth * max_depth)))
    return temperature 

    

####
# stepwise functions
####
def PlotSlabForces(filein, fileout, **kwargs):
    '''
    Plot slab surface profile
    Inputs:
        filein (str): path to input
        fileout (str): path to output
        kwargs (dict):
    '''
    assert(os.path.isfile(filein))
    ## load data: forces
    data = np.loadtxt(filein)
    depths = data[:, 0]
    buoyancies = data[:, 1]
    total_buoyancy = LA.norm(buoyancies, 1)  # total buoyancy
    buoyancie_gradients = data[:, 2]
    pressure_lower = data[:, 3]
    pressure_upper = data[:, 4]
    differiential_pressure = data[:, 5]
    differiential_pressure_v = data[:, 6]
    differiential_pressure_v_1 = data[:, 7]
    compensation = data[:, 8]
    dynamic_pressure_lower = data[:, 9]
    dynamic_pressure_upper = data[:, 10]
    differiential_dynamic_pressure = data[:, 11]
    differiential_dynamic_pressure_v = data[:, 12]
    v_zeros = np.zeros(data.shape[0])
    fig = plt.figure(tight_layout=True, figsize=(15, 15))
    gs = gridspec.GridSpec(3, 3) 
    # figure 1: show forces
    ax = fig.add_subplot(gs[0, 0]) 
    ax.plot(buoyancie_gradients, depths/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(pressure_upper, depths/1e3, 'c--', label='sigma_n (upper) (N/m2)')
    ax.plot(pressure_lower, depths/1e3, 'm--', label='sigma_n (lower) (N/m2)')
    ax.set_title('Buoyancy gradients and total pressure')
    ax.set_xlabel('Pressure (Pa)')
    ax.legend()
    ax.invert_yaxis()
    # figure 2: buoyancy and vertical pressure differences
    ax = fig.add_subplot(gs[0, 1]) 
    ax.plot(buoyancie_gradients, depths/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(differiential_pressure, depths/1e3, 'm--', label='Pressure differences (N/m2)')
    ax.plot(differiential_pressure_v, depths/1e3, 'r--', label='Vertical pressure differences (N/m2)')
    ax.plot(v_zeros, depths/1e3, 'k--')
    ax.invert_yaxis()
    ax.set_title('Buoyancy gradients and pressure differences')
    ax.set_xlabel('Pressure (Pa)')
    ax.set_ylabel('Depth (km)')
    ax.legend()
    # figure 3: buoyancy and vertical pressure differences: difined otherwise
    ax = fig.add_subplot(gs[0, 2]) 
    ax.plot(buoyancie_gradients, depths/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(differiential_pressure_v_1, depths/1e3, '--', color=mcolors.CSS4_COLORS['lightcoral'], label='Vertical pressure differences 1 (N/m2)')
    ax.plot(v_zeros, depths/1e3, 'k--')
    ax.invert_yaxis()
    ax.set_title('Buoyancy gradients and pressure differences (defined otherwise)')
    ax.set_xlabel('Pressure (Pa)')
    ax.set_ylabel('Depth (km)')
    # figure 4: field of compensation
    ax = fig.add_subplot(gs[1, 0]) 
    ax.plot(compensation, depths/1e3, 'k')
    ax.invert_yaxis()
    ax.set_title("Field of compensation")
    ax.set_xlim([-10, 10])
    ax.set_xlabel('Compensation')
    # figure 5: buoyancy and vertical pressure differences - with a 500 km limit
    mask = depths < 500e3
    ax = fig.add_subplot(gs[1, 1]) 
    ax.plot(buoyancie_gradients[mask], depths[mask]/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(differiential_pressure[mask], depths[mask]/1e3, 'm--', label='Pressure differences (N/m2)')
    ax.plot(differiential_pressure_v[mask], depths[mask]/1e3, 'r--', label='Vertical pressure differences (N/m2)')
    ax.plot(v_zeros[mask], depths[mask]/1e3, 'k--')
    ax.set_ylim([0, 500]) # set y limit
    ax.invert_yaxis()
    ax.set_title('Buoyancy gradients and pressure differences, depth in [0, 500] km')
    ax.set_xlabel('Pressure (Pa)')
    ax.set_ylabel('Depth (km)')
    ax.legend()
    # figure 6: dynamic pressure
    ax = fig.add_subplot(gs[1, 2]) 
    ax.plot(dynamic_pressure_upper, depths/1e3, 'c--', label='Dynamic P (upper) (N/m2)')
    ax.plot(dynamic_pressure_lower, depths/1e3, 'm--', label='Dynamic P (lower) (N/m2)')
    ax.set_title('Buoyancy gradients and dynamic pressure')
    ax.set_xlabel('Pressure (Pa)')
    ax.legend()
    ax.invert_yaxis()
    # figure 7: dynamic pressure differences
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(buoyancie_gradients, depths/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(differiential_dynamic_pressure, depths/1e3, 'm--', label='Dynamic P differences (N/m2)')
    ax.plot(differiential_dynamic_pressure_v, depths/1e3, 'r--', label='Vertical dynamic P differences (N/m2)')
    ax.plot(v_zeros, depths/1e3, 'k--')
    ax.set_title('Buoyancy gradients and differential dynamic pressure')
    ax.set_xlabel('Pressure (Pa)')
    ax.legend()
    ax.invert_yaxis() 
    # figure 8: dynamic pressure, in the upper 400 km
    mask = depths < 400e3
    ax = fig.add_subplot(gs[2, 1]) 
    ax.plot(buoyancie_gradients[mask], depths[mask]/1e3, 'b', label='Buoyancy gradients (N/m2)')
    ax.plot(dynamic_pressure_upper[mask], depths[mask]/1e3, 'c--', label='Dynamic P (upper) (N/m2)')
    ax.plot(dynamic_pressure_lower[mask], depths[mask]/1e3, 'm--', label='Dynamic P (lower) (N/m2)')
    ax.plot(differiential_dynamic_pressure[mask], depths[mask]/1e3, 'r--', label='Dynamic P differences (N/m2)')
    ax.plot(v_zeros[mask], depths[mask]/1e3, 'k--')
    ax.set_ylim([0, 400]) # set y limit
    ax.set_title('Buoyancy gradients and dynamic pressure')
    ax.set_xlabel('Pressure (Pa)')
    ax.legend()
    ax.invert_yaxis()
    fig.suptitle('Buoyancy (total %.4e N/m2)' % total_buoyancy)
    fig.tight_layout()
    plt.savefig(fileout)
    print("PlotSlabForces: plot figure", fileout)

def SlabMorphology(case_dir, vtu_snapshot, **kwargs):
    '''
    Wrapper for using PVTK class to get slab morphology
    Inputs:
        case_dir (str): case directory
        vtu_snapshot (int): index of file in vtu outputs
        kwargs:
            project_velocity - whether the velocity is projected to the tangential direction
    '''
    indent = kwargs.get("indent", 0)  # indentation for outputs
    findmdd = kwargs.get("findmdd", False)
    mdd_dx0 = kwargs.get('mdd_dx0', 10e3)
    mdd_dx1 = kwargs.get('mdd_dx1', 10e3)
    findmdd_tolerance = kwargs.get("findmdd_tolerance", 0.05)
    project_velocity = kwargs.get('project_velocity', False)
    mdd = -1.0 # an initial value
    print("%s%s: Start" % (indent*" ", Utilities.func_name()))
    output_slab = kwargs.get('output_slab', False)
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
    if geometry == "chunk":
        Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    else:
        Xmax = Visit_Options.options['XMAX']
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    if findmdd:
        try:
            mdd = VtkP.FindMDD(tolerance=findmdd_tolerance, dx0=mdd_dx0, dx1=mdd_dx1)
        except ValueError:
            mdd = - 1.0
    # output slab profile
    if output_slab:
        slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
        slab_internal = VtkP.ExportSlabInternal(output_xy=True)
        o_slab_env0 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env0_%05d.vtp" % (vtu_step)) # envelop 0
        o_slab_env1 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env1_%05d.vtp" % (vtu_step)) # envelop 1
        o_slab_in = os.path.join(case_dir,\
            "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
        VtkPp.ExportPolyDataFromRaw(slab_envelop0[:, 0], slab_envelop0[:, 1], None, None, o_slab_env0) # write the polydata
        # np.savetxt(o_slab_env0, slab_envelop0)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env0))
        VtkPp.ExportPolyDataFromRaw(slab_envelop1[:, 0], slab_envelop1[:, 1], None, None, o_slab_env1) # write the polydata
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env1))
        np.savetxt(o_slab_in, slab_internal)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
    # process trench, slab depth, dip angle
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    if project_velocity:
        vsp_magnitude, vov_magnitude = VtkP.ExportVelocity(project_velocity=True)
    else:
        vsp, vov = VtkP.ExportVelocity()
        vsp_magnitude = np.linalg.norm(vsp, 2)
        vov_magnitude = np.linalg.norm(vov, 2)
    # generate outputs
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e%-14.4e%-14.4e%-14.4e"\
    % (vtu_step, step, _time, trench, slab_depth, dip_100, vsp_magnitude, vov_magnitude)
    if findmdd:
        outputs += "%-14.4e" % (mdd)
    outputs += "\n"
    print("%s%s" % (indent*" ", outputs)) # debug
    return vtu_step, outputs

def SlabMorphology_dual_mdd(case_dir, vtu_snapshot, **kwargs):
    '''
    Wrapper for using PVTK class to get slab morphology, uses two distinct mdd_dx1 value
    to get both the partial and the final coupling point.
    Inputs:
        case_dir (str): case directory
        vtu_snapshot (int): index of file in vtu outputs
        kwargs:
            project_velocity - whether the velocity is projected to the tangential direction
    '''
    indent = kwargs.get("indent", 0)  # indentation for outputs
    findmdd = kwargs.get("findmdd", False)
    output_ov_ath_profile = kwargs.get("output_ov_ath_profile", False)
    output_path = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    findmdd_tolerance = kwargs.get("findmdd_tolerance", 0.05)
    depth_distant_lookup = kwargs.get("depth_distant_lookup", 200e3)
    dip_angle_depth_lookup = kwargs.get("dip_angle_depth_lookup", None)
    dip_angle_depth_lookup_interval = kwargs.get("dip_angle_depth_lookup_interval", 60e3)
    project_velocity = kwargs.get('project_velocity', False)
    find_shallow_trench = kwargs.get("find_shallow_trench", False)
    mdd = -1.0 # an initial value
    print("%s%s: Start" % (indent*" ", Utilities.func_name()))
    output_slab = kwargs.get('output_slab', None)
    # todo_crust
    n_crust = kwargs.get("n_crust", 1)
    if n_crust == 1:
        crust_fields = ['spcrust']
    elif n_crust == 2:
        crust_fields = ['spcrust_up', 'spcrust_low']
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("SlabMorphology_dual_mdd: processing %s" % filein)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    # vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    if geometry == "chunk":
        Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    else:
        Xmax = Visit_Options.options['XMAX']
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spharz', 'velocity', 'viscosity'] + crust_fields
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(crust_fields + ['spharz'], prepare_slab_distant_properties=True, depth_distant_lookup=depth_distant_lookup)
    # todo_shallow
    if find_shallow_trench:
        if n_crust == 2:
            raise NotImplementedError()
        try:
            # This is not a stable algorithm yet.
            outputs = VtkP.PrepareSlabShallow(Visit_Options.options['THETA_REF_TRENCH'])
        except Exception:
            trench_shallow = np.finfo(np.float16).tiny
        else:
            x, y, z = outputs["corrected"]["points"]
            _, _, trench_shallow = Utilities.cart2sph(x, y, z)
    if findmdd:
        try:
            mdd1 = VtkP.FindMDD(tolerance=findmdd_tolerance, dx1=-Visit_Options.options["INITIAL_SHEAR_ZONE_THICKNESS"])
            mdd2 = VtkP.FindMDD(tolerance=findmdd_tolerance, dx1=10e3)
        except ValueError:
            mdd1 = - 1.0
            mdd2 = - 1.0
    # output slab profile
    if output_slab is not None:
        slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
        slab_internal = VtkP.ExportSlabInternal(output_xy=True)
        if output_slab == "vtp":
            o_slab_env0 = os.path.join(case_dir,\
                "vtk_outputs", "slab_env0_%05d.vtp" % (vtu_step)) # envelop 0
            o_slab_env1 = os.path.join(case_dir,\
                "vtk_outputs", "slab_env1_%05d.vtp" % (vtu_step)) # envelop 1
            o_slab_in = os.path.join(case_dir,\
                "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
            VtkPp.ExportPolyDataFromRaw(slab_envelop0[:, 0], slab_envelop0[:, 1], None, None, o_slab_env0) # write the polydata
            # np.savetxt(o_slab_env0, slab_envelop0)
            print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env0))
            VtkPp.ExportPolyDataFromRaw(slab_envelop1[:, 0], slab_envelop1[:, 1], None, None, o_slab_env1) # write the polydata
            print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env1))
            np.savetxt(o_slab_in, slab_internal)
            print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
        # todo_o_env
        if output_slab == "txt":
            o_slab_env =  os.path.join(case_dir, \
                "vtk_outputs", "slab_env_%05d.txt" % (vtu_step)) # envelop 1
            slab_env_outputs = np.concatenate([slab_envelop0, slab_envelop1], axis=1) 
            slab_env_output_header = "X0 Y0 X1 Y1"
            np.savetxt(o_slab_env, slab_env_outputs, header=slab_env_output_header)
            print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env))
            o_slab_in = os.path.join(case_dir,\
                "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
            np.savetxt(o_slab_in, slab_internal)
            print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
    # output a profile distant to the slab
    n_distant_profiel_sample = 100
    v_distant_profile = None; query_viscs = None
    if output_ov_ath_profile:
        header = "# 1: x (m)\n# 2: y (m)\n# 3: depth (m)\n# 4: velocity_h (m/s)\n# 5: velocity_r (m/s)\n# 6: viscosity (pa*s)\n# 7: fixing (by 0.1 deg)\n"
        coords, depths, v_distant_profile, query_viscs, value_fixing = VtkP.ExportOvAthenProfile(depth_distant_lookup, n_sample=n_distant_profiel_sample)
        outputs = np.concatenate((coords, depths.reshape((-1, 1)), v_distant_profile,\
            query_viscs.reshape((-1, 1)), value_fixing.reshape((-1, 1))), axis=1)
        o_file = os.path.join(output_path, "ov_ath_profile_%05d.txt" % (vtu_step))
        with open(o_file, 'w') as fout:
            fout.write(header)  # output header
        with open(o_file, 'a') as fout:
            np.savetxt(fout, outputs, fmt="%20.8e")  # output data
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_file))
    # process trench, slab depth, dip angle
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    if project_velocity:
        vsp_magnitude, vov_magnitude = VtkP.ExportVelocity(project_velocity=True)
    else:
        vsp, vov = VtkP.ExportVelocity()
        vsp_magnitude = np.linalg.norm(vsp, 2)
        vov_magnitude = np.linalg.norm(vov, 2)
    # generate outputs
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e%-14.4e%-14.4e%-14.4e"\
    % (vtu_step, step, _time, trench, slab_depth, dip_100, vsp_magnitude, vov_magnitude)
    if findmdd:
        outputs += "%-14.4e %-14.4e" % (mdd1, mdd2)
    if output_ov_ath_profile:
        outputs += "%-14.4e %-14.4e" % (v_distant_profile[n_distant_profiel_sample-1, 0], query_viscs[n_distant_profiel_sample-1])
    if dip_angle_depth_lookup is not None:
        outputs += "%-14.4e" % (VtkP.GetDipAtDepth(dip_angle_depth_lookup, dip_angle_depth_lookup_interval))
    if find_shallow_trench:
        outputs += "%-14.4e" % (trench_shallow)
    
    outputs += "\n"
    print("%s%s" % (indent*" ", outputs)) # debug
    return vtu_step, outputs


def SlabAnalysis(case_dir, vtu_snapshot, o_file, **kwargs):
    '''
    Perform analysis on the slab, this would output a file including the
    buoyancy forces of the slab and the pressures on the slab surface.
    Inputs:
        kwargs(dict):
            output_slab - output slab file
            use_dT - use temperature difference as the criteria for the slab surface.
    '''
    # read in parameters
    indent = kwargs.get("indent", 0)  # indentation for outputs
    print("%s%s: Start" % (indent*" ", Utilities.func_name()))
    output_slab = kwargs.get('output_slab', False)
    output_poly_data = kwargs.get('output_poly_data', True)
    use_dT = kwargs.get('use_dT', False)
    dT = kwargs.get('dT', -100.0)
    slab_envelop_interval = kwargs.get("slab_envelop_interval", 5e3)
    ha_file = os.path.join(case_dir, "output", "depth_average.txt")
    assert(os.path.isfile(ha_file))
    output_path = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    # look for the input file named as "output/solution/solution-%05d.pvtu"
    filein = os.path.join(case_dir, "output", "solution",\
         "solution-%05d.pvtu" % (vtu_snapshot))
    assert(os.path.isfile(filein))
    # get parameters
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    # initiate class
    VtkP = VTKP(ha_file=ha_file, time=_time, slab_envelop_interval=slab_envelop_interval)
    VtkP.ReadFile(filein)
    # fields to load
    field_names = ['T', 'p', 'density', 'spcrust', 'spharz']
    has_dynamic_pressure = int(Visit_Options.options['HAS_DYNAMIC_PRESSURE']) 
    if has_dynamic_pressure == 1:
        field_names += ['nonadiabatic_pressure']
    VtkP.ConstructPolyData(field_names, include_cell_center=True, construct_Tdiff=True)
    # include a v_profile
    r0_range = [6371e3 - 2890e3, 6371e3]
    Ro = 6371e3
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    # output poly data, debug
    if output_poly_data:
        file_out = os.path.join(output_path, "processed-%05d.vtp" % vtu_snapshot)
        VtkPp.ExportPolyData(VtkP.i_poly_data, file_out)
        file_out_1 = os.path.join(output_path, "processed_center-%05d.vtp" % vtu_snapshot)
        VtkPp.ExportPolyData(VtkP.c_poly_data, file_out_1)
    # slab envelop
    if use_dT:
        VtkP.PrepareSlabByDT(slab_threshold=dT)  # slab: differential temperature
    else:
        VtkP.PrepareSlab(['spcrust', 'spharz'])  # slab: composition
    # output slab profile
    if output_slab:
        slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
        slab_internal = VtkP.ExportSlabInternal(output_xy=True)
        o_slab_env0 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env0_%05d.vtp" % (vtu_step)) # envelop 0
        o_slab_env1 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env1_%05d.vtp" % (vtu_step)) # envelop 1
        o_slab_in = os.path.join(case_dir,\
            "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
        VtkPp.ExportPolyDataFromRaw(slab_envelop0[:, 0], slab_envelop0[:, 1], None, None, o_slab_env0) # write the polydata
        VtkPp.ExportPolyDataFromRaw(slab_envelop1[:, 0], slab_envelop1[:, 1], None, None, o_slab_env1) # write the polydata
        np.savetxt(o_slab_in, slab_internal)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
    # buoyancy
    total_buoyancy, b_profile = VtkP.SlabBuoyancy(v_profile, 5e3)  # test 5e3, 50e3
    depths_o = b_profile[:, 0]  # use these depths to generate outputs
    buoyancies = b_profile[:, 1]
    buoyancy_gradients = buoyancies / (depths_o[1] - depths_o[0])  # gradient of buoyancy
    # pressure 
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()  # raw data on the envelop and output
    fileout = os.path.join(output_path, 'slab_pressures0_%05d.txt' % (vtu_step))
    depths0, thetas0, ps0= SlabPressures(VtkP, slab_envelop0, fileout=fileout, indent=4, has_dynamic_pressure=has_dynamic_pressure)  # depth, dip angle and pressure
    fileout = os.path.join(output_path, 'slab_pressures1_%05d.txt' % (vtu_step))
    depths1, thetas1, ps1 = SlabPressures(VtkP, slab_envelop1, fileout=fileout, indent=4, has_dynamic_pressure=has_dynamic_pressure)
    ps0_o = np.interp(depths_o, depths0, ps0[:, 0])  # interpolation to uniform interval
    thetas0_o = np.interp(depths_o, depths0, thetas0)
    ps0_d_o = np.interp(depths_o, depths0, ps0[:, 3])  # dynamic pressure
    ps1_o = np.interp(depths_o, depths1, ps1[:, 0])  # interpolation to uniform interval
    thetas1_o = np.interp(depths_o, depths1, thetas1)
    ps1_d_o = np.interp(depths_o, depths1, ps1[:, 3])  # dynamic pressure
    ps_o = ps0_o - ps1_o  # this has to be minus: sides of pressure are differnent on top or below.
    ps_d_o = ps0_d_o - ps1_d_o  # dynamic pressure difference
    # pvs_o = ps0_o * np.cos(thetas0_o)  - ps1_o * np.cos(thetas0_o)
    # pvs_o1 = ps0_o * np.cos(thetas0_o)  - ps1_o * np.cos(thetas1_o)  # here we cannot multiply thetas1_o, otherwise it will be zagged
    pvs_o = ps0_o / np.tan(thetas0_o)  - ps1_o / np.tan(thetas0_o)   # Right now, I am convinced this is the right way.
    pvs_d_o = ps0_d_o / np.tan(thetas0_o)  - ps1_d_o / np.tan(thetas0_o)   # vertical component of dynamic pressure differences
    pvs_o1 = ps0_o / np.tan(thetas0_o)  - ps1_o / np.tan(thetas1_o)  # here we cannot multiply thetas1_o, otherwise it will be zagged
    compensation = pvs_o / (-buoyancy_gradients)
    outputs = np.concatenate((b_profile, buoyancy_gradients.reshape((-1, 1)),\
    ps0_o.reshape((-1, 1)), ps1_o.reshape((-1, 1)),\
    ps_o.reshape((-1, 1)), pvs_o.reshape((-1, 1)), pvs_o1.reshape((-1, 1)),\
    compensation.reshape((-1, 1)), ps0_d_o.reshape((-1, 1)), ps1_d_o.reshape((-1, 1)),\
    ps_d_o.reshape((-1, 1)), pvs_d_o.reshape((-1, 1))), axis=1)
    # output data
    # all this data are outputed just to toy with the plot of buoyancy and pressure
    header = "# 1: depth (m)\n# 2: buoyancy (N/m)\n\
# 3: buoyancy gradient (Pa)\n# 4: pressure upper (Pa) \n# 5: pressure lower (Pa)\n\
# 6: differiential pressure (Pa)\n# 7: vertical differiential pressure\n\
# 8: vertical differiential pressure 1\n# 9: compensation\n\
# 10: dynamic pressure upper (Pa)\n# 11: dynamic pressure lower (Pa)\n\
# 12: differential dynamic pressure (Pa)\n# 13: vertical differential dynamic pressure (Pa)\n"
    with open(o_file, 'w') as fout:
        fout.write(header)  # output header
    with open(o_file, 'a') as fout:
        np.savetxt(fout, outputs, fmt="%20.8e")  # output data
    print("%s: write file %s" % (Utilities.func_name(), o_file))


def SlabPressures(VtkP, slab_envelop, **kwargs):
    '''
    extract slab pressures, interpolated results onto regular grid is outputed,
    original data is returned
    Inputs:
        VtkP: VTKP class
        slab_envelop: slab envelop coordinates (x and y)
    returns:
        depths: depth of points
        thetas: dip angles
        ps: pressures
    '''
    Ro = 6371e3
    fileout = kwargs.get('fileout', None)
    indent = kwargs.get('indent', 0)
    has_dynamic_pressure = kwargs.get('has_dynamic_pressure', 0)
    rs_n = 5 # resample interval
    ip_interval = 1e3  # interval for interpolation
    # resample the original envelop dataset
    n_point = slab_envelop.shape[0]
    rs_idx = range(0, n_point, rs_n)
    slab_envelop_rs = slab_envelop[np.ix_(rs_idx, [0, 1])] # use np.ix to resample, check out numpy documentation
    slab_env_polydata = VtkPp.InterpolateGrid(VtkP.i_poly_data, slab_envelop_rs, quiet=True) # note here VtkPp is module shilofue/VtkPp, while the VtkP is the class
    temp_vtk_array = slab_env_polydata.GetPointData().GetArray('p')
    env_ps  = vtk_to_numpy(temp_vtk_array)
    # dynamic pressure is outputed -> read in
    # dynamic pressure is not outputed -> use pressure - static_pressure as an estimation
    if has_dynamic_pressure == 1:
        temp_vtk_array = slab_env_polydata.GetPointData().GetArray('nonadiabatic_pressure')
        env_dps  = vtk_to_numpy(temp_vtk_array)
        print("Read in the dynamic pressures")
    else:
        env_dps = None
    # import data onto selected points
    depths = np.zeros(slab_envelop_rs.shape[0]) # data on envelop0
    ps = np.zeros((slab_envelop_rs.shape[0], 4)) # pressure, horizontal & vertical components, dynamic pressure
    thetas = np.zeros((slab_envelop_rs.shape[0], 1))
    is_first = True
    for i in range(0, slab_envelop_rs.shape[0]):
        x = slab_envelop_rs[i, 0]
        y = slab_envelop_rs[i, 1]
        theta_xy = np.arctan2(y, x)
        r = (x*x + y*y)**0.5
        depth = Ro - r  # depth of this point
        p = env_ps[i]  # pressure of this point
        p_static = VtkP.StaticPressure([r, VtkP.Ro], theta_xy, 2000)
        if env_dps is not None:
            p_d = env_dps[i]
        else:
            p_d = p - p_static # dynamic pressure, read in or compute
        depths[i] = depth
        d1 = 0.0  # coordinate differences
        d2 = 0.0
        # here we first get a dip angle
        if is_first:
            xnext = slab_envelop_rs[i+1, 0]  # coordinates of this and the last point
            ynext = slab_envelop_rs[i+1, 1]
            theta = get_dip(x, y, xnext, ynext, VtkP.geometry)
            is_first = False
        else: 
            xlast = slab_envelop_rs[i-1, 0]  # coordinates of this and the last point
            ylast = slab_envelop_rs[i-1, 1]
            theta = get_dip(xlast, ylast, x, y, VtkP.geometry) 
        thetas[i, 0] = theta
        # then we project the pressure into vertical and horizontal
        p_v = p * np.cos(theta)
        p_h = p * np.sin(theta)
        ps[i, 0] = p
        ps[i, 1] = p_h
        ps[i, 2] = p_v
        ps[i, 3] = p_d
    temp = np.concatenate((slab_envelop_rs, thetas), axis=1)
    data_env0 = np.concatenate((temp, ps), axis=1)  # assemble all the data
    c_out = data_env0.shape[1]
    # interpolate data to regular grid & prepare outputs
    start = np.ceil(depths[0]/ip_interval) * ip_interval
    end = np.floor(depths[-1]/ip_interval) * ip_interval
    n_out = int((end-start) / ip_interval)
    data_env0_out = np.zeros((n_out, c_out+1))
    depths_out = np.arange(start, end, ip_interval)
    data_env0_out[:, 0] = depths_out
    for j in range(c_out):
        data_env0_out[:, j+1] = np.interp(depths_out, depths, data_env0[:, j]) # interpolation
        header = "# 1: depth (m)\n# 2: x (m)\n# 3: y (m)\n# 4: theta_v \n# 5: p (Pa)\n# 6: p_h (Pa) \n# 7: p_v (Pa)\n"
    with open(fileout, 'w') as fout:
        fout.write(header)
    with open(fileout, 'a') as fout: 
        np.savetxt(fout, data_env0_out, fmt='%.4e\t')
    print("%s%s: write output %s" % (' '*indent, Utilities.func_name(), fileout))
    return depths, thetas[:, 0], ps  # return depths and pressures

class CmbExtractionIndexError(Exception):
    pass

# todo_temp
def SlabTemperature_old(case_dir, vtu_snapshot, ofile=None, **kwargs):
    '''
    Perform analysis on the slab, this would output a file including the
    buoyancy forces of the slab and the pressures on the slab surface.
    Inputs:
        case_dir(str) - directory of the case
        vtu_snapshot - snapshots (step plus the level of the initial adaptive refinements)
        ofile - output file of the slab temperature profile. If this is None, then no outputs are generated
        kwargs(dict):
            output_slab - output slab file
            use_dT - use temperature difference as the criteria for the slab surface.
    '''
    # assert something
    debug = kwargs.get('debug', False)
    indent = kwargs.get("indent", 0)  # indentation for outputs
    print("%s%s: Start" % (indent*" ", Utilities.func_name()))
    output_slab = kwargs.get('output_slab', False)
    output_poly_data = kwargs.get('output_poly_data', True)
    slab_envelop_interval = kwargs.get("slab_envelop_interval", 5e3)
    max_depth = kwargs.get("max_depth", 660e3)

    fix_shallow = kwargs.get("fix_shallow", False)
    ofile_surface = kwargs.get("ofile_surface", None)
    ofile_moho = kwargs.get("ofile_moho", None)
    output_path = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution",\
         "solution-%05d.pvtu" % (vtu_snapshot))
    assert(os.path.isfile(filein))
    # get parameters
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    # initiate class
    VtkP = VTKP(time=_time, slab_envelop_interval=slab_envelop_interval, slab_shallow_cutoff=25e3)
    VtkP.ReadFile(filein)
    # fields to load
    field_names = ['T', 'p', 'density', 'spcrust', 'spharz']
    has_dynamic_pressure = int(Visit_Options.options['HAS_DYNAMIC_PRESSURE']) 
    if has_dynamic_pressure == 1:
        field_names += ['nonadiabatic_pressure']
    VtkP.ConstructPolyData(field_names, include_cell_center=True, construct_Tdiff=False)
   
    # include a v_profile
    r0_range = [6371e3 - 2890e3, 6371e3]
    Ro = 6371e3
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    
    # output poly data, debug
    if output_poly_data:
        file_out = os.path.join(output_path, "processed-%05d.vtp" % vtu_snapshot)
        VtkPp.ExportPolyData(VtkP.i_poly_data, file_out)
        file_out_1 = os.path.join(output_path, "processed_center-%05d.vtp" % vtu_snapshot)
        VtkPp.ExportPolyData(VtkP.c_poly_data, file_out_1)
    
    # slab envelop
    VtkP.PrepareSlab(['spcrust', 'spharz'], prepare_moho='spcrust')  # slab: composition
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
    moho_envelop = VtkP.ExportSlabmohoCoord()
    if output_slab:
        slab_internal = VtkP.ExportSlabInternal(output_xy=True)
        o_slab_env0 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env0_%05d.vtp" % (vtu_step)) # envelop 0
        o_slab_env1 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env1_%05d.vtp" % (vtu_step)) # envelop 1
        o_moho_env = os.path.join(case_dir,\
            "vtk_outputs", "moho_env_%05d.vtp" % (vtu_step)) # envelop 0
        o_slab_in = os.path.join(case_dir,\
            "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
        VtkPp.ExportPolyDataFromRaw(slab_envelop0[:, 0], slab_envelop0[:, 1], None, None, o_slab_env0) # write the polydata
        VtkPp.ExportPolyDataFromRaw(slab_envelop1[:, 0], slab_envelop1[:, 1], None, None, o_slab_env1) # write the polydata
        # export the envelop of the core-mantle boundary
        VtkPp.ExportPolyDataFromRaw(moho_envelop[:, 0], moho_envelop[:, 1], None, None, o_moho_env) # write the polydata
        np.savetxt(o_slab_in, slab_internal)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
    rs_n = 5 # resample interval
    ip_interval = 1e3  # interval for interpolation

    # resample the origin slab surface
    n_point = slab_envelop1.shape[0]
    rs_idx = range(0, n_point, rs_n)
    slab_envelop_rs_raw = slab_envelop1[np.ix_(rs_idx, [0, 1])]  # for slab surface
    
    if fix_shallow:
        # append the shallow trench point
        outputs_shallow = VtkP.PrepareSlabShallow(Visit_Options.options['THETA_REF_TRENCH'])
        slab_envelop_rs_raw = np.vstack((outputs_shallow["corrected"]["points"][0:2], slab_envelop_rs_raw)) 
    else:
        outputs_shallow = None
    
    depths_raw = Ro - (slab_envelop_rs_raw[:, 0]**2.0 + slab_envelop_rs_raw[:, 1]**2.0)**0.5
    id_max = np.where(depths_raw < max_depth)[0][-1]
    
    depths = depths_raw[0: id_max+1]
    slab_envelop_rs = slab_envelop_rs_raw[0: id_max+1, :]

    # resample the original moho
    if debug:
        print("moho_envelop: ")  # screen outputs
        print(moho_envelop)
    try:
        moho_envelop_rs_raw = moho_envelop[np.ix_(rs_idx, [0, 1])]
        if debug:
            print("moho_envelop_rs_raw: ")  # screen outputs
            print(moho_envelop_rs_raw)
    except IndexError as e:
        rs_idx_last = rs_idx[-1]
        moho_envelop_rs_raw_length = moho_envelop.shape[0]
        raise IndexError("the last index to extract is %d, while the shape of the moho_envelop is %d" % (rs_idx_last, moho_envelop_rs_raw_length))
    
    if fix_shallow:
        # append the shallow trench point
        r, theta, phi = Utilities.cart2sph(*outputs_shallow["corrected"]["points"])
        r1 = r - Visit_Options.options["INITIAL_SHEAR_ZONE_THICKNESS"]
        moho_shallow = np.array([r1 * np.sin(theta) * np.cos(phi), r1 * np.sin(theta) * np.sin(phi), 0])
        moho_envelop_rs_raw = np.vstack((moho_shallow[0:2], moho_envelop_rs_raw))
    
    
    depths_moho_raw = Ro - (moho_envelop_rs_raw[:, 0]**2.0 + moho_envelop_rs_raw[:, 1]**2.0)**0.5
    
    depths_moho = depths_moho_raw[0: id_max+1] # match the array for the surface
    moho_envelop_rs = moho_envelop_rs_raw[0: id_max+1, :]
    
    id_valid = np.where(depths_moho > 0.0)[0][-1] # reason: max depth could be shallower than the slab surface
    
    # interpolate the curve
    # start = np.ceil(depths[0]/ip_interval) * ip_interval
    start = np.ceil(depths[0]/ip_interval) * ip_interval
    end = np.floor(depths[-1]/ip_interval) * ip_interval
    n_out = int((end-start) / ip_interval)
    depths_out = np.arange(start, end, ip_interval)

    slab_Xs = interp1d(depths, slab_envelop_rs[:, 0], kind='cubic')(depths_out)
    slab_Ys = interp1d(depths, slab_envelop_rs[:, 1], kind='cubic')(depths_out)

    mask_moho = ((depths_out > depths_moho[0]) & (depths_out < depths_moho[id_valid]))
    moho_Xs = np.zeros(depths_out.shape)
    moho_Ys = np.zeros(depths_out.shape)
    moho_Xs[mask_moho] = interp1d(depths_moho[0: id_valid+1], moho_envelop_rs[0: id_valid+1, 0], kind='cubic')(depths_out[mask_moho])
    moho_Ys[mask_moho] = interp1d(depths_moho[0: id_valid+1], moho_envelop_rs[0: id_valid+1, 1], kind='cubic')(depths_out[mask_moho])

    # interpolate the T 
    slab_env_polydata = VtkPp.InterpolateGrid(VtkP.i_poly_data, np.column_stack((slab_Xs, slab_Ys)), quiet=True) # note here VtkPp is module shilofue/VtkPp, while the VtkP is the class
    env_Ttops  = vtk_to_numpy(slab_env_polydata.GetPointData().GetArray('T'))
    
    moho_env_polydata = VtkPp.InterpolateGrid(VtkP.i_poly_data, np.column_stack((moho_Xs, moho_Ys)), quiet=True) # note here VtkPp is module shilofue/VtkPp, while the VtkP is the class
    env_Tbots = vtk_to_numpy(moho_env_polydata.GetPointData().GetArray('T'))
    
    mask = (env_Tbots < 1.0) # fix the non-sense values
    env_Tbots[mask] = -np.finfo(np.float32).max
    if debug:
        print("env_Tbots")  # screen outputs
        print(env_Tbots)

    # output 
    if ofile is not None:
        # write output if a valid path is given
        data_env0 = np.zeros((depths_out.size, 7)) # output: x, y, Tbot, Ttop
        data_env0[:, 0] = depths_out
        data_env0[:, 1] = slab_Xs
        data_env0[:, 2] = slab_Ys
        data_env0[:, 3] = moho_Xs
        data_env0[:, 4] = moho_Ys
        data_env0[:, 5] = env_Tbots
        data_env0[:, 6] = env_Ttops
        c_out = data_env0.shape[1]
        # interpolate data to regular grid & prepare outputs
        for j in range(c_out):
            header = "# 1: depth (m)\n# 2: x (m)\n# 3: y (m)\n# 4: x bot (m)\n# 5: y bot (m)\n# 6: Tbot (K)\n# 7: Ttop (K)\n"
        with open(ofile, 'w') as fout:
            fout.write(header)
        with open(ofile, 'a') as fout: 
            np.savetxt(fout, data_env0, fmt='%.4e\t')
        print("%s%s: write output %s" % (' '*indent, Utilities.func_name(), ofile))
    if ofile_surface is not None:
        # write output of surface and moho profile
        VtkPp.ExportPolyDataFromRaw(slab_Xs, slab_Ys, None, None, ofile_surface) # write the polydata
    if ofile_moho is not None:
        VtkPp.ExportPolyDataFromRaw(moho_Xs, moho_Ys, None, None, ofile_moho) # write the polydata
        
    return depths, env_Ttops, env_Tbots  # return depths and pressures


def SlabTemperature1(case_dir, vtu_snapshot, **kwargs):
    '''
    a wrapper for the SlabTemperature function
    '''
    output_path = os.path.join(case_dir, "vtk_outputs", "temperature")
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))

    if not os.path.isdir(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    ofile = os.path.join(output_path, "slab_temperature_%05d.txt" % (vtu_step))
    ofile_surface = os.path.join(output_path, "slab_surface_%05d.vtp" % (vtu_step))
    ofile_moho = os.path.join(output_path, "slab_moho_%05d.vtp" % (vtu_step))
    
    kwargs["ofile_surface"] = ofile_surface
    kwargs["ofile_moho"] = ofile_moho

    try: 
        SlabTemperature(case_dir, vtu_snapshot, ofile=ofile, **kwargs)
    except ValueError:
        warnings.warn("Generation of file for vtu_snapshot %d" % vtu_snapshot)


def WedgeT(case_dir, vtu_snapshot, **kwargs):
    '''
    export mantle temperature
    '''
    output_path = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution",\
         "solution-%05d.pvtu" % (vtu_snapshot))
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # test 1 output slab grid & envelop
    fileout = os.path.join(output_path, 'wedge_T100_%05d.txt' % (vtu_snapshot))
    VtkP.ExportWedgeT(fileout=fileout)
    assert(os.path.isfile(fileout))

def TrenchT(case_dir, vtu_snapshot, **kwargs):
    '''
    Export the trench temperature profiles
    '''
    filein = os.path.join(case_dir, "output", "solution",\
         "solution-%05d.pvtu" % (vtu_snapshot))
    output_path = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    fileout =  os.path.join(output_path, 'trench_T_%05d.txt' % (vtu_step))
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    VtkP.ExportTrenchT(fileout=fileout)
    assert(os.path.isfile(fileout))


def ShearZoneGeometry(case_dir, vtu_snapshot, **kwargs):
    indent = kwargs.get("indent", 0)  # indentation for outputs
    assert(os.path.isdir(case_dir))
    # fix the output directory
    vtk_o_dir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(vtk_o_dir):
        os.mkdir(vtk_o_dir)
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    img_o_dir = os.path.join(img_dir, "shear_zone")
    if not os.path.isdir(img_o_dir):
        os.mkdir(img_o_dir)
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("%sSlabMorphology: processing %s" % (indent*" ", filein))
    # prepare the slab
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    Dsz =  Visit_Options.options['INITIAL_SHEAR_ZONE_THICKNESS']
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # call the functions for the shear zone
    fileout = os.path.join(vtk_o_dir, "sz_%05d.txt" % vtu_step)
    VtkP.PrepareSZ(fileout, Dsz=Dsz)
    assert(os.path.isfile(fileout))  # assert file generation
    # plot
    fig_path = os.path.join(img_o_dir, "sz_thickness_%05d.png" % vtu_step) 
    fig, ax = plt.subplots()
    MorphPlotter = SLABPLOT("plot_slab")
    MorphPlotter.PlotShearZoneThickness(case_dir, axis=ax, filein=fileout, label='shear zone thickness')
    ax.legend()
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))  # assert figure generation
    print("%s%s: figure generated %s" % (indent*" ", Utilities.func_name(), fig_path))


def PlotSlabTemperature(case_dir, vtu_snapshot, **kwargs):
    '''
    Process slab envelops and plot the slab temperature
    '''
    indent = kwargs.get("indent", 0)  # indentation for outputs
    assert(os.path.isdir(case_dir))
    # read some options
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    # fix the output directory
    vtk_o_dir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(vtk_o_dir):
        os.mkdir(vtk_o_dir)
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    img_o_dir = os.path.join(img_dir, "slab_temperature")
    if not os.path.isdir(img_o_dir):
        os.mkdir(img_o_dir)
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("%sSlabMorphology: processing %s" % (indent*" ", filein))
    o_file = os.path.join(vtk_o_dir, "slab_temperature")
    if os.path.isfile(o_file):
        os.remove(o_file)
    assert(os.path.isdir(case_dir))
    _, _, _ = SlabTemperature_old(case_dir, vtu_snapshot, o_file, output_slab=True)
    assert(os.path.isfile(o_file))  # assert the outputs of slab and cmb envelops
    # plot
    fig_path = os.path.join(img_o_dir, "slab_temperature_%05d.png" % vtu_step) 
    fig, ax = plt.subplots()
    MorphPlotter = SLABPLOT("plot_slab")
    MorphPlotter.PlotSlabT(case_dir, axis=ax, filein=o_file, label='slab temperature', xlims=[273.0, 1273.0], ylims=[25e3, 250e3])
    ax.legend()
    ax.invert_yaxis()
    fig.savefig(fig_path)
    # assert figure generation and screen outputs
    assert(os.path.isfile(fig_path))  # assert figure generation
    print("%s%s: figure generated %s" % (indent*" ", Utilities.func_name(), fig_path))


def PlotSlabTemperatureCase(case_dir, **kwargs):
    '''
    Plot the slab temperature for the case
    '''
    indent = kwargs.get("indent", 0)  # indentation for outputs
    time_range = kwargs.get("time_range", None)
    plot_eclogite = kwargs.get("plot_eclogite", False)
    assert(os.path.isdir(case_dir))
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    img_o_dir = os.path.join(img_dir, "slab_temperature")
    if not os.path.isdir(img_o_dir):
        os.mkdir(img_o_dir)
    # plot
    fig, ax = plt.subplots()
    MorphPlotter = SLABPLOT("plot_slab")
    # todo_eclogite
    if plot_eclogite:
        MorphPlotter.PlotEclogite(axis=ax)
    MorphPlotter.PlotSlabTCase(case_dir, axis=ax, label='slab temperature', xlims=[273.0, 1673.0], ylims=[25e3, 250e3], time_range=time_range)
    ax.legend()
    ax.invert_yaxis()
    if time_range is None:
        fig_path = os.path.join(img_o_dir, "slab_temperature.png") 
    else:
        fig_path = os.path.join(img_o_dir, "slab_temperature_%.4eMa_%.4eMa.png" % (time_range[0]/1e6, time_range[1]/1e6)) 
    fig.savefig(fig_path)
    # assert figure generation and screen outputs
    assert(os.path.isfile(fig_path))  # assert figure generation
    print("%s%s: figure generated %s" % (indent*" ", Utilities.func_name(), fig_path))


def PlotWedgeTCase(case_dir, **kwargs):
    '''
    Plot the figure of the mantle wedge temperature
    kwargs:
        time_interval - the interval of time between two steps
    '''
    time_interval = kwargs.get("time_interval")
    ofile = os.path.join(case_dir, 'img', "wedge_T_100.png")
    SlabPlot = SLABPLOT('wedge_T')
    fig, ax = plt.subplots(figsize=(10, 4)) 
    ax, h = SlabPlot.PlotTWedge(case_dir, time_interval=time_interval, axis=ax)
    fig.colorbar(h, ax=ax, label='T (K)') 
    fig.savefig(ofile)
    assert(os.path.isfile(ofile))
    print("%s: output figure %s" % (Utilities.func_name(), ofile))

# todo_T
####
# Case-wise functions
####
# def SlabTemperatureCase(case_dir, **kwargs):
#     '''
#     run vtk and get outputs for every snapshots
#     Inputs:
#         kwargs:
#             time_interval: the interval between two processing steps
#     '''
#     # get all available snapshots
#     # the interval is choosen so there is no high frequency noises
#     time_interval_for_slab_morphology = kwargs.get("time_interval", 0.5e6)
#     offsets = kwargs.get("offsets", [])
#     Visit_Options = VISIT_OPTIONS(case_dir)
#     Visit_Options.Interpret()
#     # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
#     available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
#     print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
#     # get where previous session ends
#     vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
#     if not os.path.isdir(vtk_output_dir):
#         os.mkdir(vtk_output_dir)
#     # Initiation Wrapper class for parallel computation
#     ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_temperature', SlabTemperature1, if_rewrite=True, assemble=False, output_poly_data=False, fix_shallow=True, offsets=offsets)
#     ParallelWrapper.configure(case_dir)  # assign case directory
#     # Remove previous file
#     print("%s: Delete old slab_temperature.txt file." % Utilities.func_name())
#     ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite
#     num_cores = multiprocessing.cpu_count()
#     # loop for all the steps to plot
#     Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_snapshot)\
#     for pvtu_snapshot in available_pvtu_snapshots)  # first run in parallel and get stepwise output
#     ParallelWrapper.clear()
#     # for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
#     #    ParallelWrapper(pvtu_snapshot)
#     # pvtu_steps_o, outputs = ParallelWrapper.assemble()

def SlabMorphologyCase(case_dir, **kwargs):
    '''
    run vtk and get outputs for every snapshots
    Inputs:
        kwargs:
            rewrite: if rewrite previous results
            project_velocity - whether the velocity is projected to the tangential direction
            file_tag - apply a tag to file name, default is false
    '''
    # todo_o_env
    findmdd = kwargs.get('findmdd', False)
    findmdd_tolerance = kwargs.get('findmdd_tolerance', 0.05)
    project_velocity = kwargs.get('project_velocity', False)
    # todo_shallow
    find_shallow_trench = kwargs.get("find_shallow_trench", False)
    # todo_parallel
    use_parallel = kwargs.get('use_parallel', False)
    file_tag = kwargs.get('file_tag', False)
    output_ov_ath_profile = kwargs.get('output_ov_ath_profile', False)
    kwargs["if_rewrite"] = True
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
    # file name
    if file_tag == 'interval' and abs(time_interval_for_slab_morphology - 5e5)/5e5 > 1e-6:
        slab_morph_file_name = 'slab_morph_t%.2e.txt' % time_interval_for_slab_morphology
    else:
        slab_morph_file_name = 'slab_morph.txt'
    slab_morph_file = os.path.join(vtk_output_dir, slab_morph_file_name)
    # Initiation Wrapper class for parallel computation
    # ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_morph', SlabMorphology_dual_mdd, if_rewrite=True, findmdd=findmdd, project_velocity=project_velocity, findmdd_tolerance=findmdd_tolerance)
    # parse the number of crust from the Visit_Options variable
    kwargs['n_crust'] = Visit_Options.options["N_CRUST"]
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_morph', SlabMorphology_dual_mdd, **kwargs)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file
    if os.path.isfile(slab_morph_file):
        print("%s: Delete old slab_morph.txt file." % Utilities.func_name())
        os.remove(slab_morph_file)  # delete slab morph file
    ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite
    ParallelWrapper.set_pvtu_steps(available_pvtu_snapshots)
    num_cores = multiprocessing.cpu_count()
    # loop for all the steps to plot, the parallel version doesn't work for now
    if use_parallel:
        # raise NotImplementedError("Parallel for the function %s is not properly implemented yet" % Utilities.func_name())
        Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_snapshot)\
        for pvtu_snapshot in available_pvtu_snapshots)  # first run in parallel and get stepwise output
        print("call assemble_parallel")  # debug
        pvtu_steps_o, outputs = ParallelWrapper.assemble_parallel()
    else:
        for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
            ParallelWrapper(pvtu_snapshot)
        pvtu_steps_o, outputs = ParallelWrapper.assemble()
    ParallelWrapper.clear()
    # last, output
    # header
    file_header = "# 1: pvtu_step\n# 2: step\n# 3: time (yr)\n# 4: trench (rad)\n# 5: slab depth (m)\n\
# 6: 100km dip (rad)\n# 7: subducting plate velocity (m/yr)\n# 8: overiding plate velocity (m/yr)\n"
    n_col = 8
    if findmdd:
        file_header += "# 9: mechanical decoupling depth1 (m)\n# 10: mechanical decoupling depth2 (m)\n"
        n_col += 2
    if output_ov_ath_profile:
        file_header += "# %d: athenosphere velocity (m/yr)\n# %d: athenosphere viscosity (pa*s)\n" % (n_col+1, n_col+2)
        n_col += 2
    if find_shallow_trench:
        file_header += "# %d: shallow trench (rad)\n" % (n_col+1)
        n_col += 1
    # output file name
    appendix = ""
    if abs(time_interval_for_slab_morphology - 0.5e6) / 0.5e6 > 1e-6:
        appendix = "_t%.2e" % time_interval_for_slab_morphology
    output_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph' + appendix + '.txt')
    # output data
    if not os.path.isfile(output_file):
        with open(output_file, 'w') as fout:
            fout.write(file_header)
            for output in outputs:
                fout.write("%s" % output)
        print('Created output: %s' % output_file)
    else:
        with open(output_file, 'a') as fout:
            for output in outputs:
                fout.write("%s" % output)
        print('Updated output: %s' % output_file)


def WedgeTCase(case_dir, **kwargs):
    '''
    run vtk and get outputs for every snapshots
    Inputs:
        kwargs:
            time_interval: the interval between two processing steps
    '''
    # get all available snapshots
    # the interval is choosen so there is no high frequency noises
    time_interval = kwargs.get("time_interval", 0.5e6)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval)
    print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
    # get where previous session ends
    vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(vtk_output_dir):
        os.mkdir(vtk_output_dir)
    # Initiation Wrapper class for parallel computation
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('wedgeT', WedgeT, if_rewrite=True, assemble=False, output_poly_data=False)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file
    print("%s: Delete old slab_temperature.txt file." % Utilities.func_name())
    ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite
    num_cores = multiprocessing.cpu_count()
    # loop for all the steps to plot
    Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_snapshot)\
    for pvtu_snapshot in available_pvtu_snapshots)  # first run in parallel and get stepwise output
    ParallelWrapper.clear()
    # for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
    #    ParallelWrapper(pvtu_snapshot)
    # pvtu_steps_o, outputs = ParallelWrapper.assemble()


def TrenchTCase(case_dir, **kwargs):
    '''
    run vtk and get outputs of trench temperature profile for every snapshots
    Inputs:
        case_dir: the directory of case
        kwargs:
            time_interval: the interval between two processing steps
    '''
    # get all available snapshots
    # the interval is choosen so there is no high frequency noises
    time_interval = kwargs.get("time_interval", 0.5e6)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval)
    print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
    # get where previous session ends
    vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(vtk_output_dir):
        os.mkdir(vtk_output_dir)
    # Initiation Wrapper class for parallel computation
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('trenchT', TrenchT, if_rewrite=True, assemble=False, output_poly_data=False)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file
    print("%s: Delete old slab_temperature.txt file." % Utilities.func_name())
    ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite
    num_cores = multiprocessing.cpu_count()
    # loop for all the steps to plot
    Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_snapshot)\
    for pvtu_snapshot in available_pvtu_snapshots)  # first run in parallel and get stepwise output
    ParallelWrapper.clear()
    # for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
    #    ParallelWrapper(pvtu_snapshot)
    # pvtu_steps_o, outputs = ParallelWrapper.assemble()
    

def PlotSlabForcesCase(case_dir, vtu_step, **kwargs):
    '''
    Inputs:
        case_dir (str): case directory
        step : step to plot
        kwargs(dict):
            output_slab - output slab file
    '''
    output_slab = kwargs.get('output_slab', False)
    assert(os.path.isdir(case_dir))
    Visit_Options = VISIT_OPTIONS(case_dir)
    # call function
    Visit_Options.Interpret()
    vtu_snapshot = int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) + vtu_step
    vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(vtk_output_dir):
        os.mkdir(vtk_output_dir)
    ofile = os.path.join(vtk_output_dir, "slab_forces_%05d" % vtu_step)
    SlabAnalysis(case_dir, vtu_snapshot, ofile, output_slab=output_slab)
    # plot figure
    img_dir = os.path.join(case_dir, 'img')
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    fig_ofile = os.path.join(img_dir, "slab_forces_%05d.png" % vtu_step)
    PlotSlabForces(ofile, fig_ofile)



def PlotSlabShape(in_dir, vtu_step):
    '''
    Plot the shape of the slab, debug usage
    Inputs:
        in_dir (str): directory containing the data file
            a. a "slab_env0_{vtu_step}.txt" and a "slab_env1_{vtu_step}.txt" file
            b. a "slab_internal_{vtu_step}.txt" file
        vtu_step: step in visualization.
    '''
    fig, ax = plt.subplots()
    file_env0 = os.path.join(in_dir, "slab_env0_%05d.txt" % vtu_step)
    file_env1 = os.path.join(in_dir, "slab_env1_%05d.txt" % vtu_step)
    file_inter = os.path.join(in_dir, "slab_internal_%05d.txt" % vtu_step)
    slab_env0 = np.loadtxt(file_env0)
    slab_env1 = np.loadtxt(file_env1)
    slab_inter = np.loadtxt(file_inter)
    ax.plot(slab_env0[:, 0], slab_env0[:, 1], 'b.')
    ax.plot(slab_env1[:, 0], slab_env1[:, 1], 'c.')
    ax.plot(slab_inter[:, 0], slab_inter[:, 1], 'r.')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_aspect('equal', adjustable='box')
    plt.show()


def ShearZoneGeometryCase(case_dir, **kwargs):
    '''
    plot the shear zone geometry for a single case
    '''
    indent = kwargs.get("indent", 0)  # indentation for outputs
    time_start = kwargs.get("time_start", 0.0)
    time_interval = kwargs.get("time_interval", 5e6)
    time_end = kwargs.get("time_end", 60e6)
    assert(os.path.isdir(case_dir))
    # fix the output directory
    vtk_o_dir = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(vtk_o_dir):
        os.mkdir(vtk_o_dir)
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    img_o_dir = os.path.join(img_dir, "shear_zone")
    if not os.path.isdir(img_o_dir):
        os.mkdir(img_o_dir)
    # initialize the VTK object
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_start=time_start, time_interval=time_interval, time_end=time_end)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    Dsz =  Visit_Options.options['INITIAL_SHEAR_ZONE_THICKNESS']
    # initiate the function and object for the figure
    fig, ax = plt.subplots()
    MorphPlotter = SLABPLOT("plot_slab")
    length = len(available_pvtu_snapshots)
    normalizer = [ float(i)/(length-1) for i in range(length) ] 
    colors = cm.rainbow(normalizer) 
    i = 0
    for vtu_snapshot in available_pvtu_snapshots:
        # prepare the slab
        vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
        _time, step = Visit_Options.get_time_and_step(vtu_step)
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
        if not os.path.isfile(filein):
            raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
        else:
            print("%sSlabMorphology: processing %s" % (indent*" ", filein))
        VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
        VtkP.ReadFile(filein)
        field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
        VtkP.ConstructPolyData(field_names, include_cell_center=True)
        VtkP.PrepareSlab(['spcrust', 'spharz'])
        # call the functions for the shear zone
        fileout = os.path.join(vtk_o_dir, "sz_%05d.txt" % vtu_step)
        VtkP.PrepareSZ(fileout, Dsz=Dsz)
        assert(os.path.isfile(fileout))  # assert file generation
        plot_initial = False
        if i == 0:
            plot_initial = True
        MorphPlotter.PlotShearZoneThickness(case_dir, plot_initial, axis=ax, filein=fileout,\
                                            label="t = %.4e Ma" % (_time/1e6), color=colors[i])
        i += 1
    # plot
    ax.legend()
    fig_path = os.path.join(img_o_dir, "sz_thickness_combined_s%.1f_i%.1f_e%.1f.png"\
        % (time_start/1e6, time_interval/1e6, time_end/1e6)) 
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))  # assert figure generation
    print("%s%s: figure generated %s" % (indent*" ", Utilities.func_name(), fig_path))


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

    class MorphFileReadingError(Exception):
        pass
 
    def ReadWedgeT(self, case_dir, **kwargs):
        '''
        read the wedge_T100 files and rearrange data
        '''
        time_interval = kwargs.get("time_interval", 0.5e6)
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        initial_adaptive_refinement = int(self.prm['Mesh refinement']['Initial adaptive refinement'])
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        elif geometry == 'box':
            Ro = float(self.prm['Geometry model']['Box']['Y extent'])
        else:
            raise ValueError('Invalid geometry')
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval)
        i = 0
        for pvtu_step in available_pvtu_snapshots:
            file_in_path = os.path.join(case_dir, 'vtk_outputs', 'wedge_T100_%05d.txt' % pvtu_step)
            print("file_in_path: ", file_in_path)  # debug
            Utilities.my_assert(os.access(file_in_path, os.R_OK), FileExistsError, "File %s doesn\'t exist" % file_in_path)
            self.ReadHeader(file_in_path)
            self.ReadData(file_in_path)
            col_x = self.header['x']['col']
            col_y = self.header['y']['col']
            col_T = self.header['T']['col']
            xs = self.data[:, col_x]
            ys = self.data[:, col_y]
            if i == 0: 
                rs = (xs**2.0 + ys**2.0)**0.5
                depthes = Ro - rs # compute depth
                # Ts = np.zeros((depthes.size, max_pvtu_step - min_pvtu_step + 1))
                Ts = np.zeros((depthes.size, len(available_pvtu_snapshots)))
            Ts[:, i] = self.data[:, col_T]
            i += 1
        return depthes, Ts

    # todo_shallow
    def PlotMorph(self, case_dir, **kwargs):
        '''
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        save_pdf = kwargs.get("save_pdf", False) 
        compare_shallow_trench = kwargs.get("compare_shallow_trench", False)
        # initialization
        findmdd = False
        # path
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_dir = os.path.join(img_dir, 'morphology')
        if not os.path.isdir(morph_dir):
            os.mkdir(morph_dir)
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        time_interval = times[1] - times[0]

        try: 
            col_pvtu_shallow_trench = self.header['shallow_trench']['col']
        except KeyError:
            col_pvtu_shallow_trench = None
            shallow_trenches = None
        else:
            shallow_trenches = self.data[:, col_pvtu_shallow_trench]

        if time_interval < 0.5e6:
            warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
            if shallow_trenches is not None:
                shallow_trenches_migration_length = (shallow_trenches - shallow_trenches[0]) * Ro  # length of migration
            else:
                shallow_trenches_migration_length = None
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
            shallow_trenches_migration_length = None # not implemented yet
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        if shallow_trenches_migration_length is not None:
            shallow_trench_velocities = np.gradient(shallow_trenches_migration_length, times)
        else:
            shallow_trench_velocities = None
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]
        try:
            col_mdd1 = self.header['mechanical_decoupling_depth1']['col']
            col_mdd2 = self.header['mechanical_decoupling_depth2']['col']
        except KeyError:
            pass
        else:
            findmdd = True
            mdds1 = self.data[:, col_mdd1]
            mdds2 = self.data[:, col_mdd2]
        # trench velocity
        # start figure
        if findmdd:
            _size = (20, 10)
            gs = gridspec.GridSpec(4, 2) 
        else:
            _size = (15, 10)
            gs = gridspec.GridSpec(3, 2) 
        fig = plt.figure(tight_layout=True, figsize=(15, 10)) 
        fig.subplots_adjust(hspace=0)
        # 1: trench & slab movement
        ax = fig.add_subplot(gs[0, 0:2]) 
        ax_tx = ax.twinx()
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color='tab:orange', label='trench position (km)')
        if shallow_trenches_migration_length is not None:
            lns0_1 = ax.plot(times/1e6, shallow_trenches_migration_length/1e3, '--', color='tab:orange')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylabel('Trench Position (km)', color="tab:orange")
        ax.tick_params(axis='x', labelbottom=False) # labels along the bottom edge are off
        ax.tick_params(axis='y', labelcolor="tab:orange")
        ax.grid()
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, 'k-', label='slab depth (km)')
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        # ax.legend(lns, labs)
        # 2: velocity
        ax = fig.add_subplot(gs[1, 0:2]) 
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color='tab:orange', label='trench velocity (cm/yr)')
        if shallow_trenches_migration_length is not None:
            lns0_1 = ax.plot(times/1e6, shallow_trench_velocities*1e2, '--', color='tab:orange')
        lns1 = ax.plot(times/1e6, sp_velocities*1e2, '-', color='tab:blue', label='subducting plate (cm/yr)')
        lns2 = ax.plot(times/1e6, ov_velocities*1e2, '-', color='tab:purple', label='overiding velocity (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='sinking velocity (cm/yr)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylim((-10, 10))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.set_xlabel('Times (Myr)')
        ax.grid()
        ax.legend()
        # 2.1: velocity smaller, no y limit, to show the whole curve
        ax = fig.add_subplot(gs[2, 0:2]) 
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color="tab:orange", label='trench velocity (cm/yr)')
        lns1 = ax.plot(times/1e6, sp_velocities*1e2, '-', color='tab:blue', label='subducting plate (cm/yr)')
        lns2 = ax.plot(times/1e6, ov_velocities*1e2, '-', color='tab:purple', label='overiding velocity (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='trench velocity (cm/yr)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylabel('Velocity (whole, cm/yr)')
        ax.grid()
        # 3: the mechanical decoupling depth, only if the data is found
        if findmdd:
            ax = fig.add_subplot(gs[3, 0:2]) 
            ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
            lns3 = ax.plot(times/1e6, mdds1/1e3, '-', color="tab:blue", label='mdd1 (km)')
            lns4 = ax.plot(times/1e6, mdds2/1e3, '-', color="c", label='mdd2 (km)')
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
            ax.set_ylabel('Depth (km)')
            ax.grid()
            ax.legend()
        fig.tight_layout()
        # save figure
        if save_pdf:
            o_path = os.path.join(morph_dir, 'trench.pdf')
            plt.savefig(o_path, format="pdf", bbox_inches="tight")
        else:
            o_path = os.path.join(morph_dir, 'trench.png')
            plt.savefig(o_path)
        print("%s: figure %s generated" % (Utilities.func_name(), o_path))
        if compare_shallow_trench:
            # compare shallow trench to trench
            assert(shallow_trenches_migration_length is not None)
            # Create a figure with 3 subplots arranged vertically
            fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 8))
            # plot location
            lns0 = ax0.plot(times/1e6, trenches, '-', color='tab:orange', label='trench')
            lns0_1 = ax0.plot(times/1e6, shallow_trenches, '--', color='tab:orange', label="shallow trench")
            ax0.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
            ax0.set_ylabel('Trench Position')
            ax0.grid()
            ax0.legend()
            lns1 = ax1.plot(times/1e6, trench_velocities*1e2, '-', color='tab:orange', label='trench velocity (cm/yr)')
            lns1_1 = ax1.plot(times/1e6, shallow_trench_velocities*1e2, '--', color='tab:orange')
            ax1.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
            ax1.set_ylabel('Velocity (whole, cm/yr)')
            ax1.set_xlabel('Times (Myr)')
            ax1.grid()
            o_path = os.path.join(morph_dir, 'shallow_trench_compare.pdf')
            plt.savefig(o_path)
            print("%s: figure %s generated" % (Utilities.func_name(), o_path))
    
    def PlotMorphAnime(self, case_dir, **kwargs):
        '''
        Plot slab morphology for animation
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            time -a time to plot
        '''
        _time = kwargs.get('time', '0.0')
        # path
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_dir = os.path.join(img_dir, 'morphology')
        if not os.path.isdir(morph_dir):
            os.mkdir(morph_dir)
        temp_dir = os.path.join(morph_dir, 'temp')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        time_interval = times[1] - times[0]
        if time_interval < 0.5e6:
            warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]

        # 1: trench & slab movement
        gs = gridspec.GridSpec(2, 1) 
        fig = plt.figure(tight_layout=True, figsize=(10, 10)) 
        ax = fig.add_subplot(gs[0, 0]) 
        ax_tx = ax.twinx()
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color='tab:orange', label='trench position (km)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ylim0 = np.floor(np.min(trenches_migration_length)/1e5)*1e2
        ylim1 = np.floor(np.max(trenches_migration_length)/1e5)*1e2
        ax.set_ylim((ylim0, ylim1))
        ax.set_ylabel('Trench Position (km)', color="tab:orange")
        ax.set_xlabel('Time (Ma)')
        ax.tick_params(axis='y', labelcolor="tab:orange")
        ax.grid()
        temp_ts = _time * np.ones(100)
        temp_ys = np.linspace(ylim0, ylim1, 100)
        ax.plot(temp_ts/1e6, temp_ys, 'c--') # plot a vertical line
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, 'k-', label='slab depth (km)')
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        ax.legend(lns, labs, loc='lower right')
        # ax1: velocity
        ax = fig.add_subplot(gs[1, 0]) 
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color="tab:orange", label='trench velocity (cm/yr)')
        lns1 = ax.plot(times/1e6, sp_velocities*1e2, '-', color='tab:blue', label='subducting plate (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='trench velocity (cm/yr)')
        temp_ts = _time * np.ones(100)
        temp_ys = np.linspace(-10, 10, 100)
        ax.plot(temp_ts/1e6, temp_ys, 'c--') # plot a vertical line
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylim((-10, 10))
        ax.xaxis.set_ticks_position("top")
        ax.set_ylabel('Velocity (whole, cm/yr)')
        ax.grid()
        ax.legend(loc='lower right')
        # save figure
        o_path = os.path.join(temp_dir, 'trench_t%.4e.png' % _time)
        fig.savefig(o_path)
        print("%s: save figure %s" % (Utilities.func_name(), o_path))
        plt.close()

    def PlotMorphPublication(self, case_dir, **kwargs):
        '''
        Plot slab morphology for publication
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            time -a time to plot
        '''
        time_interval = kwargs.get('time_interval', 5e6)
        time_range = kwargs.get('time_range', None)
        if time_range is not None:
            assert(len(time_range) == 2)
        time_markers = kwargs.get("time_markers", [])
        vlim = kwargs.get("vlim", [-10, 10])
        save_pdf = kwargs.get("save_pdf", False)
        assert(len(vlim) == 2)
        # path
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_dir = os.path.join(img_dir, 'morphology')
        if not os.path.isdir(morph_dir):
            os.mkdir(morph_dir)
        # visit option class 
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        # input file name
        appendix = ""
        if abs(time_interval - 0.5e6) / 0.5e6 > 1e-6:
            appendix = "_t%.2e" % time_interval
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph' + appendix + '.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        col_athenosphere_velocity = self.header['athenosphere_velocity']['col']
        col_athenosphere_viscosity = self.header['athenosphere_viscosity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        time_interval = times[1] - times[0]
        if time_interval < 0.5e6:
            warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]
        athenosphere_velocities = self.data[:, col_athenosphere_velocity]
        athenosphere_viscosities = self.data[:, col_athenosphere_viscosity]

        # 1: trench & slab movement
        gs = gridspec.GridSpec(4, 1) 
        fig = plt.figure(tight_layout=True, figsize=(20, 40)) 
        ax = fig.add_subplot(gs[0, 0]) 
        ax_tx = ax.twinx()
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color='tab:orange', label='trench position (km)')
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ylim0 = np.floor(np.min(trenches_migration_length)/1e5)*1e2
        ylim1 = np.floor(np.max(trenches_migration_length)/1e5)*1e2
        ax.set_ylim((ylim0, ylim1))
        ax.set_ylabel('Trench Position (km)', color="tab:orange")
        ax.set_xlabel('Time (Ma)')
        ax.tick_params(axis='y', labelcolor="tab:orange")
        ax.grid()
        for _time in time_markers:
            temp_ts = _time * np.ones(100)
            temp_ys = np.linspace(ylim0, ylim1, 100)
            ax.plot(temp_ts/1e6, temp_ys, 'c--', dashes=(10, 10), alpha=0.7) # plot a vertical line
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, 'k-', label='slab depth (km)')
        if time_range is None:
            ax_tx.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax_tx.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        ax.legend(lns, labs, loc='lower right')
        # ax1: velocity
        # ax1, part 1: velcoity
        ax = fig.add_subplot(gs[1, 0]) 
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        ln_v_tr = ax.plot(times/1e6, trench_velocities*1e2, '-', color="tab:orange", label='trench velocity (cm/yr)')
        ln_v_sub = ax.plot(times/1e6, sp_velocities*1e2, '-', color='tab:blue', label='subducting plate (cm/yr)')
        ln_v_ov = ax.plot(times/1e6, ov_velocities*1e2, '-', color='tab:purple', label='overiding plate (cm/yr)')
        ln_v_ath = ax.plot(times/1e6, athenosphere_velocities*1e2, '-', color='c', label='athenosphere (cm/yr)')
        ln_v_sink = ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='sink velocity (cm/yr)')
        for _time in time_markers:
            temp_ts = _time * np.ones(200)
            temp_ys = np.linspace(-100, 100, 200)
            ax.plot(temp_ts/1e6, temp_ys, 'c--', dashes=(10, 10), alpha=0.7) # plot a vertical line
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ax.set_xlabel('Time (Ma)')
        ax.set_ylim((vlim[0], vlim[1]))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.grid()
        # ax1, part 1: trench position
        ax_tx = ax.twinx()
        ln_tr = ax_tx.plot(times/1e6, trenches_migration_length/1e3, '*', color='tab:orange', label='trench position (km)')
        ax_tx.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ylim0 = np.floor(np.min(trenches_migration_length)/1e5)*1e2
        ylim1 = np.floor(np.max(trenches_migration_length)/1e5)*1e2
        if time_range is None:
            ax_tx.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax_tx.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ax_tx.set_ylim((ylim0, ylim1))
        ax_tx.set_ylabel('Trench Position (km)', color="tab:orange")
        ax_tx.tick_params(axis='y', labelcolor="tab:orange")
        lns = ln_tr + ln_v_tr + ln_v_sub + ln_v_sink + ln_v_ov + ln_v_ath
        labs = [I.get_label() for I in lns]
        ax.legend(lns, labs, loc='upper right')
        # read athenosphere dataset
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=float(time_interval))
        depth_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        time_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        viscosity_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        velocity_h_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        for i in range(len(available_pvtu_snapshots)):
            vtu_snapshot = available_pvtu_snapshots[i]
            _time, _ = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
            vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
            slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'ov_ath_profile_%.5d.txt' % vtu_step)
            assert(os.path.isfile(slab_morph_file))
            data = np.loadtxt(slab_morph_file)
            depths = data[:, 2]
            depth_mesh[i, :] = depths
            time_mesh[i, :] = np.ones(100) * _time
            velocities_h = data[:, 3]
            velocity_h_mesh[i, :] = velocities_h
            viscosities = data[:, 5]
            viscosity_mesh[i, :] = viscosities
        # plot the viscosity
        ax = fig.add_subplot(gs[2, 0])
        h = ax.pcolormesh(time_mesh/1e6, depth_mesh/1e3, np.log10(viscosity_mesh), cmap=ccm.roma,\
            vmin=np.log10(Visit_Options.options["ETA_MIN"]), vmax=np.log10(Visit_Options.options["ETA_MAX"]))
        ax.invert_yaxis()
        fig.colorbar(h, ax=ax, label='log(viscosity (Pa*s))', orientation="horizontal")
        ax.set_xlabel("Time (Ma)")
        ax.set_ylabel("Depth (km)")
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        # plot the velocity
        ax = fig.add_subplot(gs[3, 0])
        # todo_2dmorph
        vlim_for_ath = kwargs.get("vlim_for_ath", None)
        v_min = np.floor(velocity_h_mesh.min()*100.0 / 5.0) * 5.0
        v_max = np.ceil(velocity_h_mesh.max()*100.0 / 5.0) * 5.0
        if vlim_for_ath is not None:
            assert(type(vlim_for_ath) == list and len(vlim_for_ath) == 2)
            v_min = vlim_for_ath[0]
            v_max = vlim_for_ath[1]
        h = ax.pcolormesh(time_mesh/1e6, depth_mesh/1e3, velocity_h_mesh*100.0, cmap=ccm.vik, vmin=v_min, vmax=v_max)
        ax.invert_yaxis()
        fig.colorbar(h, ax=ax, label='velocity (cm/yr)', orientation="horizontal")
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Time (Ma)")
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        # save figure
        o_path = os.path.join(morph_dir, 'trench_t%.2e' % time_interval)
        fig.savefig(o_path + '.png')
        print("%s: save figure %s" % (Utilities.func_name(), o_path + '.png'))
        if save_pdf:
            fig.savefig(o_path + '.pdf')
            print("%s: save figure %s" % (Utilities.func_name(), o_path + '.pdf'))
        plt.close()

    def PlotMorphPublicationBillen18(self, case_dir, **kwargs):
        '''
        Plot slab morphology for publication
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            time -a time to plot
        '''
        time_interval = kwargs.get('time_interval', 5e6)
        time_range = kwargs.get('time_range', None)
        time_markers = kwargs.get("time_markers", [])
        vlim = kwargs.get("vlim", [-10, 10])
        save_pdf = kwargs.get("save_pdf", False)
        assert(len(vlim) == 2)
        # path
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_dir = os.path.join(img_dir, 'morphology')
        if not os.path.isdir(morph_dir):
            os.mkdir(morph_dir)
        # visit option class 
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        # input file name
        appendix = ""
        if abs(time_interval - 0.5e6) / 0.5e6 > 1e-6:
            appendix = "_t%.2e" % time_interval
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph' + appendix + '.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        col_athenosphere_velocity = self.header['athenosphere_velocity']['col']
        col_athenosphere_viscosity = self.header['athenosphere_viscosity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        time_interval = times[1] - times[0]
        if time_interval < 0.5e6:
            warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]
        athenosphere_velocities = self.data[:, col_athenosphere_velocity]
        athenosphere_viscosities = self.data[:, col_athenosphere_viscosity]

        # start and end time
        time0, time1 = times[0]/1e6, times[-1]/1e6
        if time_range is not None:
            assert(len(time_range) == 2)
            time0, time1 = time_range[0] / 1e6, time_range[1] / 1e6

        # ax1, part 1: velcoity
        gs = gridspec.GridSpec(3, 1) 
        fig = plt.figure(tight_layout=True, figsize=(20, 30)) 
        ax = fig.add_subplot(gs[0, 0]) 
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        ln_v_tr = ax.plot(times/1e6, trench_velocities*1e2, '-', color="tab:orange", label='trench velocity (cm/yr)')
        ln_v_sub = ax.plot(times/1e6, sp_velocities*1e2, '-', color='tab:blue', label='subducting plate (cm/yr)')
        ln_v_ov = ax.plot(times/1e6, ov_velocities*1e2, '-', color='tab:purple', label='overiding plate (cm/yr)')
        ln_v_ath = ax.plot(times/1e6, athenosphere_velocities*1e2, '-', color='r', label='athenosphere (cm/yr)')
        for _time in time_markers:
            temp_ts = _time * np.ones(200)
            temp_ys = np.linspace(-100, 100, 200)
            ax.plot(temp_ts/1e6, temp_ys, 'c--', dashes=(10, 10), alpha=0.7) # plot a vertical line
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ax.set_ylim((vlim[0], vlim[1]))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.tick_params(axis='x', which='both', direction='in', labelbottom=False)
        ax.grid()
        lns = ln_v_tr + ln_v_sub + ln_v_ov + ln_v_ath
        labs = [I.get_label() for I in lns]
        ax.legend(lns, labs, loc='upper right')
        # read athenosphere dataset
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=float(time_interval))
        depth_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        time_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        viscosity_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        velocity_h_mesh = np.zeros([len(available_pvtu_snapshots), 100])
        for i in range(len(available_pvtu_snapshots)):
            vtu_snapshot = available_pvtu_snapshots[i]
            _time, _ = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
            vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
            slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'ov_ath_profile_%.5d.txt' % vtu_step)
            # Utilities.my_assert(os.path.isfile(slab_morph_file), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), slab_morph_file))
            if os.path.isfile(slab_morph_file):
                data = np.loadtxt(slab_morph_file)
                depths = data[:, 2]
                depth_mesh[i, :] = depths
                time_mesh[i, :] = np.ones(100) * _time
                velocities_h = data[:, 3]
                velocity_h_mesh[i, :] = velocities_h
                viscosities = data[:, 5]
                viscosity_mesh[i, :] = viscosities
            else:
                # a method to fix invalid timestep:
                # make sure this plots outside of the figure
                depths = np.ones(100) * (-1)
                time_mesh[i, :] = np.ones(100) * _time
                velocities_h = data[:, 3]
        # plot the viscosity
        ax = fig.add_subplot(gs[1, 0])
        h = ax.pcolormesh(time_mesh/1e6, depth_mesh/1e3, np.log10(viscosity_mesh), cmap=ccm.roma,\
            vmin=np.log10(Visit_Options.options["ETA_MIN"]), vmax=np.log10(Visit_Options.options["ETA_MAX"]))
        ax.invert_yaxis()
        fig.colorbar(h, ax=ax, label='log(viscosity (Pa*s))', orientation="horizontal")
        ax.set_ylabel("Depth (km)")
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        ax.tick_params(axis='x', which='both', direction='in', labelbottom=False)
        # plot the velocity
        ax = fig.add_subplot(gs[2, 0])
        # todo_2dmorph
        vlim_for_ath = kwargs.get("vlim_for_ath", None)
        v_min = np.floor(velocity_h_mesh.min()*100.0 / 5.0) * 5.0
        v_max = np.ceil(velocity_h_mesh.max()*100.0 / 5.0) * 5.0
        if vlim_for_ath is not None:
            assert(type(vlim_for_ath) == list and len(vlim_for_ath) == 2)
            v_min = vlim_for_ath[0]
            v_max = vlim_for_ath[1]
        h = ax.pcolormesh(time_mesh/1e6, depth_mesh/1e3, velocity_h_mesh*100.0, cmap=ccm.vik, vmin=v_min, vmax=v_max)
        ax.invert_yaxis()
        fig.colorbar(h, ax=ax, label='velocity (cm/yr)', orientation="horizontal")
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Time (Ma)")
        if time_range is None:
            ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        else:
            ax.set_xlim((time_range[0]/1e6, time_range[1]/1e6))  # set x limit
        fig.subplots_adjust(hspace=0)
        # save figure
        o_path = os.path.join(morph_dir, 'trench_b18_t%.2e' % time_interval)
        fig.savefig(o_path + '.png')
        print("%s: save figure %s" % (Utilities.func_name(), o_path + '.png'))
        if save_pdf:
            fig.savefig(o_path + '.pdf')
            print("%s: save figure %s" % (Utilities.func_name(), o_path + '.pdf'))
        plt.close()

    def PlotTWedge(self, case_dir, **kwargs):
        '''
        plot the mantle wedge temperature on top of the 100-km deep slab.
        '''
        time_interval = kwargs.get("time_interval", 0.5e6)
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        depthes, Ts = self.ReadWedgeT(case_dir)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
        times = []
        for snapshot in available_pvtu_snapshots:
            _time, _ = Visit_Options.get_time_and_step_by_snapshot(snapshot)
            times.append(_time)
        times = np.array(times)
        tt, dd = np.meshgrid(times, depthes)
        h = ax.pcolormesh(tt/1e6,dd/1e3,Ts, shading='gouraud') 
        ax.invert_yaxis()
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_xlabel('Times (Myr)')
        ax.set_ylabel('Depth (km)')
        return ax, h

    
    def PlotTrenchVelocity(self, case_dir, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        label_all = kwargs.get('label_all', False)
        color = kwargs.get('color', None)
        time_range = kwargs.get('time_range', [])
        v_range = kwargs.get('v_range', [])
        fix_v_range = kwargs.get('fix_v_range', False)
        if label_all:
            # if label_all, append labels, otherwise don't
            labels = ["trench velocity", "subducting plate", "overiding velocity", "sinking velocity"]
        else:
            labels = [None, None, None, None]
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        Utilities.my_assert(os.path.isfile(slab_morph_file), FileExistsError, "%s doesn't exist" % slab_morph_file)
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]
        # trench velocity
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color=color, label=labels[0])
        lns1 = ax.plot(times/1e6, sp_velocities*1e2, ':', color=color, label=labels[1])
        lns2 = ax.plot(times/1e6, ov_velocities*1e2, '-.', color=color, label=labels[2])
        ax.plot(times/1e6, sink_velocities*1e2, '--', color=color, label=labels[3])
        if time_range != []:
            xlims = time_range
        else:
            xlims = (np.min(times), np.max(times))
        ax.set_xlim(xlims[0]/1e6, xlims[1]/1e6)  # set x limit
        # for the limit of y, there are 3 options: a. fix_v_range would give a (-10, 20);
        # b. assigne a v_range will apply that value; c. by default, the min value of 
        # the trench velocity and the max value of the subducting velocity will be used.
        if v_range != []:
            ylims = v_range
        else:
            mask = (times > xlims[0]) & (times < xlims[1])
            ylims = [-0.15, np.max(sp_velocities[mask])]
        if fix_v_range:
            ax.set_ylim((-10, 20))
        else:
            ax.set_ylim((ylims[0]*1e2, ylims[1]*1e2))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.set_xlabel('Times (Myr)')
        ax.grid()
        ax.legend()
        # lns = lns0 + lns1
        # labs = [I.get_label() for I in lns]
        # return lns, labs
    
    
    def PlotTrenchPosition(self, case_dir, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        label = kwargs.get('label', [None, None])
        assert(len(label) == 2)
        color = kwargs.get('color', None)
        time_range = kwargs.get('time_range', [])
        tp_range = kwargs.get('tp_range', [])
        sd_range = kwargs.get('sd_range', [])
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        Utilities.my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,\
        "prm file %s cannot be opened" % prm_file)
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        Utilities.my_assert(os.path.isfile(slab_morph_file), FileExistsError, "%s doesn't exist" % slab_morph_file)
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        col_pvtu_sp_v = self.header['subducting_plate_velocity']['col']
        col_pvtu_ov_v = self.header['overiding_plate_velocity']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        sp_velocities = self.data[:, col_pvtu_sp_v]
        ov_velocities = self.data[:, col_pvtu_ov_v]
        # trench velocity
        # 1: trench & slab movement
        ax_tx = ax.twinx()
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color=color, label=label[0])
        if time_range != []:
            xlims = time_range
        else:
            xlims = (np.min(times), np.max(times))
        ax.set_xlim(xlims[0]/1e6, xlims[1]/1e6)  # set x limit
        ax.set_xlabel("Time (Myr)")
        if tp_range != []:
            ylims = tp_range
        else:
            ylims = (np.min(trenches_migration_length), np.max(trenches_migration_length))
        ax.set_ylim(ylims[0]/1e3, ylims[1]/1e3)
        ax.set_ylabel('Trench Position (km)')
        ax.tick_params(axis='x') # labels along the bottom edge are off
        ax.tick_params(axis='y')
        ax.grid()
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, '--', color=color, label=label[1])
        if sd_range != []:
            ylims = sd_range
        else:
            ylims = (np.min(slab_depthes), np.max(slab_depthes))
        ax_tx.set_ylim(ylims[0]/1e3, ylims[1]/1e3)
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        return lns, labs

    def PlotShearZoneThickness(self, case_dir, plot_initial=True, **kwargs):
        '''
        a variation of the PlotMorph function: used for combining results
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        filein = kwargs.get("filein", None)
        if filein is not None:
            sz_file = filein
        else:
            sz_file = os.path.join(case_dir, 'vtk_outputs', 'shear_zone.txt')
        Utilities.my_assert(os.path.isfile(sz_file), FileExistsError, "%s doesn't exist" % sz_file)
        label = kwargs.get('label', None)
        xlims = kwargs.get('xlims', None)
        ylims = kwargs.get('ylims', None)
        _color = kwargs.get("color", "tab:blue")
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        wb_file = os.path.join(case_dir, 'case.wb')
        assert(os.access(wb_file, os.R_OK))
        # get the initial thickness of the shear zone
        with open(wb_file, 'r') as fin:
            wb_dict = json.load(fin)
        i0 = FindWBFeatures(wb_dict, 'Subducting plate')
        sp_dict = wb_dict['features'][i0]
        initial_thickness = sp_dict["composition models"][0]["max depth"]
        # read data
        self.ReadHeader(sz_file)
        self.ReadData(sz_file)
        col_depth = self.header['depth']['col']
        col_theta_min = self.header['theta_min']['col']
        col_theta_max = self.header['theta_max']['col']
        depths = self.data[:, col_depth]
        theta_mins = self.data[:, col_theta_min]
        theta_maxs = self.data[:, col_theta_max]
        # convert to thickness along strike
        num = depths.size
        thicks = np.zeros(num)
        for i in range(num):
            r_min = Ro - depths[i]
            theta_min = theta_mins[i]
            theta_max = theta_maxs[0]
            thick = 0.0
            for j in range(0, num):
                r_max = Ro - depths[j]
                theta_max = theta_maxs[j]
                thick_temp = Utilities.point2dist([theta_min, r_min], [theta_max, r_max], geometry)
                if j == 0:
                    thick = thick_temp
                if thick_temp < thick:
                    thick = thick_temp
            thicks[i] = thick 
        # plot
        mask = (depths > initial_thickness) & (theta_mins > 0.0) # points with theta min < 0.0 are those on the surface
        ax.plot(depths[mask]/1e3, thicks[mask]/1e3, label=label, color=_color) # debug
        if plot_initial:
            ax.plot(depths/1e3, initial_thickness*np.ones(depths.size)/1e3, 'k--')
        if xlims is not None:
            # set x limit
            assert(len(xlims) == 2)
            ax.set_xlim([xlims[0]/1e3, xlims[1]/1e3])
        if ylims is not None:
            # set y limit
            assert(len(ylims) == 2)
            ax.set_ylim([ylims[0]/1e3, ylims[1]/1e3])
        ax.set_xlabel("Depth (km)")
        ax.set_ylabel("Thickness (km)")

    def PlotSlabT(self, case_dir, **kwargs):
        '''
        Plot the slab temperature
        Inputs:
            case_dir (str) - directory of case
            kwargs(dict):
                axis - a matplotlib axis
        '''
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        filein = kwargs.get("filein", None)
        if filein is not None:
            temp_file = filein
        else:
            temp_file = os.path.join(case_dir, 'vtk_outputs', 'shear_zone.txt')
        Utilities.my_assert(os.path.isfile(temp_file), FileExistsError, "%s doesn't exist" % temp_file)
        label = kwargs.get('label', None)
        xlims = kwargs.get('xlims', None)
        ylims = kwargs.get('ylims', None)
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        # read data
        self.ReadHeader(temp_file)
        self.ReadData(temp_file)
        col_depth = self.header['depth']['col']
        col_Tbot = self.header['Tbot']['col']
        col_Ttop = self.header['Ttop']['col']
        depths = self.data[:, col_depth]
        Tbots = self.data[:, col_Tbot]
        Ttops = self.data[:, col_Ttop]
        ax.plot(Ttops, depths/1e3, label=label, color="tab:blue")
        mask = (Tbots > 0.0)  # non-sense values are set to negative when these files are generated
        ax.plot(Tbots[mask], depths[mask]/1e3, label=label, color="tab:green")
        if xlims is not None:
            # set temperature limit
            assert(len(xlims) == 2)
            ax.set_xlim([xlims[0], xlims[1]])
        if ylims is not None:
            # set depth limit
            assert(len(ylims) == 2)
            ax.set_ylim([ylims[0]/1e3, ylims[1]/1e3])
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Depth (km)")

    def PlotSlabTCase(self, case_dir, **kwargs):
        '''
        Plot the temperature profile for a case by assembling results
        from individual steps
            kwargs(dict):
                axis - a matplotlib axis
                debug - print debug message
                time_range - range of the time for plotting the temperature
        '''
        # initiate
        n_plot = 100 # points in the plot
        max_plot_depth = 250e3
        ax = kwargs.get('axis', None)
        debug = kwargs.get('debug', False)
        use_degree = kwargs.get('use_degree', False)
        if ax == None:
            raise ValueError("Not implemented")
        label = kwargs.get('label', None)
        xlims = kwargs.get('xlims', None)
        ylims = kwargs.get('ylims', None)
        time_range = kwargs.get('time_range', None)
        if time_range is not None:
            assert(len(time_range) == 2)
            assert(time_range[0] < time_range[1])
        # options for slab temperature outputs
        time_interval_for_slab_morphology = kwargs.get("time_interval", 0.5e6)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
        # assert all files exist
        for pvtu_snapshot in available_pvtu_snapshots:
            vtu_step = max(0, int(pvtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
            output_path = os.path.join(case_dir, "vtk_outputs")
            temp_file = os.path.join(output_path, "slab_temperature_%05d.txt" % (vtu_step))
            assert(os.access(temp_file, os.R_OK))
        pDepths = np.linspace(0, max_plot_depth, n_plot)
        # derive the range of the slab temperatures
        pTtops_min = np.ones(n_plot) * 1e31
        pTtops_max = np.ones(n_plot) * (-1e31)
        pTtops_med = np.zeros(n_plot)
        pTtops_wt = np.zeros(n_plot)
        pTbots_min = np.ones(n_plot) * 1e31
        pTbots_max = np.ones(n_plot) * (-1e31)
        pTbots_med = np.zeros(n_plot)
        pTbots_wt = np.zeros(n_plot)
        time_last = 0.0
        for pvtu_snapshot in available_pvtu_snapshots:
            vtu_step = max(0, int(pvtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
            _time, step = Visit_Options.get_time_and_step(vtu_step)
            if time_range is not None:
                if _time < time_range[0] or _time > time_range[1]:
                    # if the step is out of the time range, skip it.
                    continue
            dtime = (_time - time_last)
            time_last = _time
            output_path = os.path.join(case_dir, "vtk_outputs")
            temp_file = os.path.join(output_path, "slab_temperature_%05d.txt" % (vtu_step))
            # read data from file
            self.ReadHeader(temp_file)
            self.ReadData(temp_file)
            col_depth = self.header['depth']['col']
            col_Tbot = self.header['Tbot']['col']
            col_Ttop = self.header['Ttop']['col']
            depths = self.data[:, col_depth]
            Tbots = self.data[:, col_Tbot]
            Ttops = self.data[:, col_Ttop]
            Tbot_func = interp1d(depths, Tbots, assume_sorted=True) 
            Ttop_func = interp1d(depths, Ttops, assume_sorted=True) 
            for i in range(n_plot):
                pDepth = pDepths[i]
                if pDepth < depths[0] or pDepth > depths[-1]:
                    # the range in the slab temperature file
                    # could be limited, so that invalid values
                    # are skipped 
                    continue
                pTtop = Ttop_func(pDepth)
                if pTtop > 0.0:
                    if pTtop < pTtops_min[i]:
                        pTtops_min[i] = pTtop
                    if pTtop > pTtops_max[i]:
                        pTtops_max[i] = pTtop
                    pTtops_med[i] += pTtop * dtime
                    pTtops_wt[i] += dtime
                pTbot = Tbot_func(pDepth)
                if pTbot > 0.0:
                    # only deal with valid values
                    if pTbot < pTbots_min[i]:
                        pTbots_min[i] = pTbot
                    if pTbot > pTbots_max[i]:
                        pTbots_max[i] = pTbot
                    pTbots_med[i] += pTbot * dtime
                    pTbots_wt[i] += dtime
        pTtops_med /= pTtops_wt
        pTbots_med /= pTbots_wt
        if debug:
            print("pTbots_min: ")  # screen outputs
            print(pTbots_min)
            print("pTbots_max: ") 
            print(pTbots_max)
            print("pTbots_med: ")
            print(pTbots_med)
        # plot result of slab surface temperature
        mask = ((pTbots_min > 0.0) & (pTbots_min < 1e4))
        ax.plot(pTbots_min[mask], pDepths[mask]/1e3, "--", label="cmb T", color="tab:green")
        mask = ((pTbots_max > 0.0) & (pTbots_max < 1e4))
        ax.plot(pTbots_max[mask], pDepths[mask]/1e3, "--", color="tab:green")
        mask = ((pTbots_med > 0.0) & (pTbots_med < 1e4))
        ax.plot(pTbots_med[mask], pDepths[mask]/1e3, "-",  color="tab:green")
        # plot result of cmb temperature
        mask = ((pTtops_min > 0.0) & (pTtops_min < 1e4))
        ax.plot(pTtops_min[mask], pDepths[mask]/1e3, "--", label="surface T", color="tab:blue")
        mask = ((pTtops_max > 0.0) & (pTtops_max < 1e4))
        ax.plot(pTtops_max[mask], pDepths[mask]/1e3, "--", color="tab:blue")
        mask = ((pTtops_med > 0.0) & (pTtops_med < 1e4))
        ax.plot(pTtops_med[mask], pDepths[mask]/1e3, "-",  color="tab:blue")
        ax.set_xlim(xlims)
        ax.set_ylim([ylims[0]/1e3, ylims[1]/1e3])
        ax.set_ylabel("Depth (km)")
        if use_degree:
            ax.set_xlabel("Temperature (C)")
        else:
            ax.set_xlabel("Temperature (K)")

    # todo_eclogite
    def PlotEclogite(self, **kwargs):
        '''
        Plot the temperature profile for a case by assembling results
        from individual steps
            kwargs(dict):
                axis - a matplotlib axis
                debug - print debug message
                time_range - range of the time for plotting the temperature
        '''
        ax = kwargs.get('axis', None)
        use_degree = kwargs.get('use_degree', False)
        p_to_depth = 3.3/100 # 100 km 3.3 GPa
        file_stern_2001 = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', 'reference', 'eclogite_stern_2001.txt')
        assert(os.path.isfile(file_stern_2001))
        file_hu_2022 = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', 'reference', 'eclogite_hernandez-uribe_2022.txt')
        assert(os.path.isfile(file_hu_2022))
        # data from stern 2001
        self.ReadHeader(file_stern_2001)
        self.ReadData(file_stern_2001)
        col_T = self.header['temperature']['col']
        col_depth = self.header['depth']['col']
        Ts = self.data[:, col_T] # C
        depths = self.data[:, col_depth]
        if not use_degree:
            Ts += 273.0 # K
        ax.plot(Ts, depths, 'k-.', label='stern_2001')
        # data from the Hernandez-uribe_2022
        self.ReadHeader(file_hu_2022)
        self.ReadData(file_hu_2022)
        col_T = self.header['temperature']['col']
        col_P = self.header['pressure']['col']
        Ts = self.data[:, col_T] # C
        depths = self.data[:, col_P] / p_to_depth
        if not use_degree:
            Ts += 273.0 # K
        ax.plot(Ts, depths, 'c-.', label='Hernandez-uribe_2022')

    def FitTrenchT(self, case_dir, vtu_snapshot):
        '''
        fit the trench temperature
        '''
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        Utilities.my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,
                'BASH_OPTIONS.__init__: case prm file - %s cannot be read' % prm_file)
        with open(prm_file, 'r') as fin:
            idict = ParseFromDealiiInput(fin)
        potential_temperature = idict.get('Adiabatic surface temperature', 1673.0)
        # read the temperature and depth profile
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
        trench_T_file = os.path.join(case_dir, 'vtk_outputs', 'trench_T_%05d.txt' % (vtu_step))
        assert(os.path.isfile(trench_T_file))
        self.ReadHeader(trench_T_file)
        self.ReadData(trench_T_file)
        col_depth = self.header['depth']['col']
        col_T = self.header['T']['col']
        depths = self.data[:, col_depth]
        Ts = self.data[:, col_T]
        tFitFunc = T_FIT_FUNC(depths, Ts, potential_temperature=1573.0)
        x0 = [1.0]  # variables: age; scaling: 40 Ma
        res = minimize(tFitFunc.PlateModel, x0)
        age_myr = res.x[0] * 40.0 # Ma
        print('step: ', vtu_step, 'age_myr: ', age_myr)
        return age_myr

    def GetSlabMorph(self, case_dir):
        '''
        read the slab_morph file
        '''
        morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        assert(os.path.isfile(morph_file))
        self.ReadHeader(morph_file)
        self.ReadData(morph_file)
        col_sp_velocity = self.header['subducting_plate_velocity']['col']
        col_ov_velocity = self.header['overiding_plate_velocity']['col']
        col_dip = self.header['100km_dip']['col']
        col_time = self.header['time']['col']
        col_trench = self.header['trench']['col']
        times = self.data[:, col_time]
        sp_velocities = self.data[:, col_sp_velocity]
        ov_velocities = self.data[:, col_ov_velocity]
        trenches = self.data[:, col_trench]
        dips = self.data[:, col_dip]
        conv_velocities = sp_velocities - ov_velocities
        return times, sp_velocities, ov_velocities, dips, trenches

    def GetAgeTrench(self, case_dir, use_thermal=True, **kwargs):
        '''
        get the ages of the subducting plate at the trench
        Inputs:
            kwargs:
                time_interval - interval between steps
        '''
        time_interval = kwargs.get('time_interval', 0.5e6)
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        if use_thermal:
            # use thermal option would fit the temperature at the trench for the individual steps
            # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
            available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval)
            print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug
            age_trenchs = []
            times = []
            for pvtu_snapshot in available_pvtu_snapshots:
                vtu_step = max(0, int(pvtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
                output_path = os.path.join(case_dir, "vtk_outputs")
                temp_file = os.path.join(output_path, "trench_T_%05d.txt" % (vtu_step))
                assert(os.access(temp_file, os.R_OK))
            for pvtu_snapshot in available_pvtu_snapshots:
                vtu_step = max(0, int(pvtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
                _time, _ = Visit_Options.get_time_and_step(vtu_step)
                times.append(_time)
                age_trench = self.FitTrenchT(case_dir, pvtu_snapshot)
                age_trenchs.append(age_trench)
            times = np.array(times)
            age_trenchs = np.array(age_trenchs)
        else:
            morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
            self.ReadHeader(morph_file)
            self.ReadData(morph_file)
            assert(os.path.isfile(morph_file))
            # read the geometry & Ro
            prm_file = os.path.join(case_dir, 'output', 'original.prm')
            assert(os.access(prm_file, os.R_OK))
            self.ReadPrm(prm_file)
            # read parameters
            geometry = self.prm['Geometry model']['Model name']
            if geometry == 'chunk':
                Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
            elif geometry == 'box':
                Ro =  float(self.prm['Geometry model']['box']["Y extent"])
            # read the spreading velocity 
            wb_file = os.path.join(case_dir, 'case.wb')
            assert(os.path.isfile(wb_file))
            with open(wb_file, 'r') as fin:
                wb_dict = json.load(fin)
            i0 = FindWBFeatures(wb_dict, 'Subducting plate')
            sp_dict = wb_dict['features'][i0]
            trench_ini = sp_dict['coordinates'][2][0] * np.pi / 180.0
            spreading_velocity = sp_dict['temperature models'][0]['spreading velocity']
            # read data in slab_morph.txt
            col_time = self.header['time']['col']
            col_trench = self.header['trench']['col']
            col_sp_velocity = self.header['subducting_plate_velocity']['col']
            times = self.data[:, col_time]
            trenchs = self.data[:, col_trench]
            sp_velocities = self.data[:, col_sp_velocity]
            # correction on the trench coordinates
            trench_correction = trench_ini - trenchs[0]
            trenchs = trenchs + trench_correction
            age_trenchs = []
            is_first = True
            dist = 0.0
            time_last = 0.0
            for i in range(times.size):
                _time = times[i]
                sp_velocity = sp_velocities[i]
                trench = trenchs[i]
                dtime = _time - time_last 
                if is_first:
                    is_first = False
                else:
                    dist += sp_velocity * dtime
                time_last = _time
                # first part is the initial age of the material at the trench
                if geometry == 'box':
                    age = (trench - dist) / spreading_velocity + _time
                elif geometry == 'chunk':
                    age = (trench * Ro - dist) / spreading_velocity + _time
                age_trenchs.append(age)
            age_trenchs = np.array(age_trenchs)
        return times, age_trenchs


    def PlotTrenchAge(self, case_dir, **kwargs):
        '''
        plot the age of the trench
        Inputs:
            kwargs:
                time_interval - interval between steps
        '''
        time_interval = kwargs.get('time_interval', 0.5e6)
        ax = kwargs.get('axis', None)
        use_thermal = kwargs.get('use_thermal', True)
        # options for slab temperature outputs
        times, age_trenchs = self.GetAgeTrench(case_dir, use_thermal, time_interval=time_interval)
        if use_thermal:
            figure_label = "age of the trench (%s)" % "thermally interpreted"
        else:
            figure_label = "age of the trench (%s)" % "motion reconstruction"
        ax.plot(times/1e6, age_trenchs, label=figure_label)
        ax.set_xlabel("Time (Ma)")
        ax.set_ylabel("Trench Age (Ma)")
    
    def PlotThermalParameter(self, case_dir, **kwargs):
        '''
        plot the age of the trench
        Inputs:
            kwargs:
                time_interval - interval between steps
                time_stable - begining of the stable subduction
                plot_velocity - plot the velocity alongside the thermal parameter
        '''
        time_interval = kwargs.get('time_interval', 0.5e6)
        time_stable = kwargs.get('time_stable', None)
        ax = kwargs.get('axis', None)
        use_thermal = kwargs.get('use_thermal', True)
        plot_velocity = kwargs.get('plot_velocity', False)
        plot_dip = kwargs.get('plot_dip', False)
        # options for slab temperature outputs
        if use_thermal:
            # not implemented yet, the array of the age_trenchs imported this way could
            # be different in size from other arrays
            raise NotImplementedError()
        _, age_trenchs = self.GetAgeTrench(case_dir, use_thermal, time_interval=time_interval)
        times, sp_velocities, ov_velocities, dips, _ = self.GetSlabMorph(case_dir)
        conv_velocities = sp_velocities - ov_velocities
        thermal_parameters = age_trenchs * conv_velocities * np.sin(dips)
        if use_thermal:
            figure_label = "thermal parameter (%s)" % "thermally interpreted"
        else:
            figure_label = "thermal parameter (%s)" % "motion reconstruction"
        ax.plot(times/1e6, thermal_parameters/1e3, label=figure_label)
        ax.set_xlabel("Time (Ma)")
        ax.set_ylabel("Thermal Parameter (km)", color='tab:blue')
        if time_stable is not None:
            # focus on the range of thermal parameter in the stable regem
            mask = (times > time_stable)
            # the ValueError would be induced if the model hasn't reached
            # the stage of stable subduction
            try:
                ymax = np.ceil(np.max(thermal_parameters[mask]) / 1e6) * 1e6
                ymin = np.floor(np.min(thermal_parameters[mask]) / 1e6) * 1e6
                ax.set_ylim([ymin/1e3, ymax/1e3])
            except ValueError:
                pass
        # converging velocity
        if plot_velocity:
            ax1 = ax.twinx()
            ax1.plot(times/1e6, conv_velocities/1e3, '--', label="conv velocity", color='tab:green')
            ax1.set_ylabel("Velocity (m/yr)", color='tab:green')
            if time_stable is not None:
                # focus on the range of thermal parameter in the stable regem
                mask = (times > time_stable)
                try:
                    ymax1 = np.ceil(np.max(conv_velocities[mask]) / 1e-3) * 1e-3
                    ymin1 = np.floor(np.min(conv_velocities[mask]) / 1e-3) * 1e-3
                    ax1.set_ylim([ymin1/1e3, ymax1/1e3])
                except ValueError:
                    pass
        if plot_dip:
            ax2 = ax.twinx()
            ax2.plot(times/1e6, np.sin(dips), '--', label="sin(dip angle)", color='tab:red')
            ax2.set_ylabel("sin(dip angle)", color='tab:red')
            if time_stable is not None:
                # focus on the range of thermal parameter in the stable regem
                mask = (times > time_stable)
                try:
                    ymax2 = np.ceil(np.max(np.sin(dips[mask])) / 1e-2) * 1e-2
                    ymin2 = np.floor(np.min(np.sin(dips[mask])) / 1e-2) * 1e-2
                    ax2.set_ylim([ymin2, ymax2])
                except ValueError:
                    pass

    def GetTimeDepthTip(self, case_dir, query_depth, **kwargs):
        '''
        todo_t660
        Get the time the slab tip is at a certain depth
        Inputs:
            case_dir (str): case directory
        '''
        filename = kwargs.get("filename", "slab_morph.txt")

        assert(os.path.isdir(case_dir))
        morph_file = os.path.join(case_dir, 'vtk_outputs', filename)
        Utilities.my_assert(os.path.isfile(morph_file), self.SlabMorphFileNotExistError, "%s is not a file." % morph_file)
        self.ReadHeader(morph_file)
        self.ReadData(morph_file)

        try: 
            col_time = self.header['time']['col']
            times = self.data[:, col_time]
            col_slab_depth = self.header['slab_depth']['col']
            slab_depths = self.data[:, col_slab_depth]
        except IndexError as e:
            # in case the file cannot be read, just return an invalid value
            return -1.0
            # raise SLABPLOT.MorphFileReadingError("Error while reading slab morphology file %s" % morph_file) from e
        query_time = -1.0
        for i in range(len(times)-1):
            _time = times[i]
            depth = slab_depths[i]
            next_depth = slab_depths[i+1]
            if depth < query_depth and next_depth > query_depth:
                next_time = times[i+1]
                query_time = (query_depth - depth) / (next_depth - depth) * next_time +\
                    (query_depth - next_depth) / (depth - next_depth) * _time
                
        return query_time

    def GetAverageVelocities(self, case_dir, t0, t1, **kwargs):
        '''
        Inputs:
            compute the average velocities between t0 and t1
        '''
        assert(os.path.isdir(case_dir))
        filename = kwargs.get("filename", "slab_morph.txt")

        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        Utilities.my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,\
        "prm file %s cannot be opened" % prm_file)
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        else:
            Ro = -1.0  # in this way, wrong is wrong
        
        # read morph file
        morph_file = os.path.join(case_dir, 'vtk_outputs', filename)
        Utilities.my_assert(os.path.isfile(morph_file), self.SlabMorphFileNotExistError, "%s is not a file." % morph_file)
        self.ReadHeader(morph_file)
        self.ReadData(morph_file)
        
        # assign initial values 
        V_sink_avg = -1.0
        V_plate_avg = -1.0
        V_ov_plate_avg = -1.0
        V_trench_avg = -1.0

        try: 
            col_time = self.header['time']['col']
            times = self.data[:, col_time]
        except IndexError:
            return V_sink_avg, V_plate_avg, V_ov_plate_avg, V_trench_avg
            
        # if the time range is invalid (t1 should be bigger), return a state of -2.0 
        if (t1 <= t0):
            V_sink_avg = -2.0
            V_plate_avg = -2.0
            V_ov_plate_avg = -2.0
            V_trench_avg = -2.0
            return V_sink_avg, V_plate_avg, V_ov_plate_avg, V_trench_avg
        
        # calculate the velocity is suitable value of t1 is provided
        if t1 < times[times.size-1]:
            col_trench = self.header['trench']['col']
            trenches = self.data[:, col_trench]
            col_slab_depth = self.header['slab_depth']['col']
            slab_depths = self.data[:, col_slab_depth]
            col_sp_velocity = self.header['subducting_plate_velocity']['col']
            col_op_velocity = self.header['overiding_plate_velocity']['col']
            sp_velocities = self.data[:, col_sp_velocity]
            op_velocities = self.data[:, col_op_velocity]
    
            # trench migration 
            trench0 = np.interp(t0, times, trenches)
            trench1 = np.interp(t1, times, trenches)
            trenches_migration_length = 0.0
            if geometry == "chunk":
                trenches_migration_length = (trench1 - trench0) * Ro  # length of migration
            elif geometry == 'box':
                trenches_migration_length = trench1 - trench0
            V_trench_avg = trenches_migration_length / (t1 - t0)
    
            # slab depth 
            slab_depth0 = np.interp(t0, times, slab_depths)
            slab_depth1 = np.interp(t1, times, slab_depths)
            V_sink_avg = (slab_depth1 - slab_depth0) / (t1 - t0)
    
            # plate motion
            # first compute the velocity of the subducting plate
            v_temp = 0.0
            w_temp = 0.0
            for i in range(times.size-1):
                time0 = times[i]
                time1 = times[i+1]
                if time0 > t0 and time1 < t1:
                    v_temp += (time1 - time0) * (sp_velocities[i] + sp_velocities[i+1])/2.0
                    w_temp += (time1 - time0)
            V_plate_avg = v_temp / w_temp
            # then compute the velocity of the overiding plate 
            v_temp = 0.0
            w_temp = 0.0
            for i in range(times.size-1):
                time0 = times[i]
                time1 = times[i+1]
                if time0 > t0 and time1 < t1:
                    v_temp += (time1 - time0) * (op_velocities[i] + op_velocities[i+1])/2.0
                    w_temp += (time1 - time0)
            V_ov_plate_avg = v_temp / w_temp

        return V_sink_avg, V_plate_avg, V_ov_plate_avg, V_trench_avg


    class SlabMorphFileNotExistError(Exception):
        pass

    def write_csv(self, case_dir, **kwargs):
        '''
        using the pandas interface to convert to csv
        Inputs:
            case_dir (str): direction of the case
        '''
        # read data
        o_csv_path = kwargs.get("o_path", None)
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        if not self.HasData():
            print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
            return 1
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        # read the geometry & Ro
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        elif geometry == 'box':
            Ro =  float(self.prm['Geometry model']['box']["Y extent"])
        # get the trench migration length 
        if geometry == "chunk":
            trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        elif geometry == 'box':
            trenches_migration_length = trenches - trenches[0]
        else:
            raise ValueError('Invalid geometry')
        # collect data 
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        # assemble in an output
        o_csv_array = np.zeros([self.data.shape[0], 2])
        o_csv_array[:, 0] = times
        o_csv_array[:, 1] = trenches
        # uses a default path in vtk_outputs if no option is giving
        if o_csv_path is None:
            o_csv_path = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.csv')
        # export
        # TODO: add field names and write R functions to parse the result
        df = pd.DataFrame(o_csv_array)
        df.to_csv(o_csv_path)


class SLABMATERIAL(LINEARPLOT): 
    '''
    A class defined to plot the Slab materials (crust and harzburgite)
    '''
    def __init__(self, _name):
        LINEARPLOT.__init__(self, _name)
        self.ha_reader =  DEPTH_AVERAGE_PLOT('DepthAverage') # reader of the horizontal average file

    def ReadFile(self, case_dir):
        '''
        Inputs:
        '''
        # read inputs
        prm_file = os.path.join(case_dir, 'output', 'original.prm')
        assert(os.access(prm_file, os.R_OK))
        self.ReadPrm(prm_file)
        # read parameters
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            self.Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
            self.geometry_extent = float(self.prm['Geometry model']['Chunk']['Chunk maximum longitude'])\
                - float(self.prm['Geometry model']['Chunk']["Chunk minimum longitude"])
            self.geometry_input = "cartesian" # this is the input that goes into the integretion function
        elif geometry == 'box':
            self.Ro =  float(self.prm['Geometry model']['box']["Y extent"])
            self.geometry_extent =  float(self.prm['Geometry model']['box']["X extent"])
            self.geometry_input = "spherical"
        else:
            return ValueError("%s: geometry must be chunk or box" % Utilities.func_name())
        # read data
        case_output_dir = os.path.join(case_dir, 'output')
        depth_average_file = os.path.join(case_output_dir, 'depth_average.txt')
        assert(os.access(depth_average_file, os.R_OK))
        self.ha_reader.ReadHeader(depth_average_file)  # inteprate header information
        self.ha_reader.ReadData(depth_average_file)  # read data
        # self.ha_reader.ManageUnits()  # mange unit to output
        self.ha_reader.SplitTimeStep()  # split time step data
   
   
    def PlotSlabMaterial(self, _time, ax):
        '''
        plot the slab material
        Inputs:
            ax: an axis is passed in for plotting
        '''
        data_list0, step = self.ha_reader.ExportDataByTime(_time, ["depth", "spcrust", "spharz"])
        depths = data_list0[:, 0]
        spcrusts = data_list0[:, 1]
        spcrust_integretions, spcrust_segmentations = self.ha_reader.GetIntegrateArray(_time, "spcrust",\
            2, self.geometry_input, self.geometry_extent, Ro=self.Ro)
        # plot
        mask = (spcrusts > 0.0)
        # ax.semilogx(spcrusts[mask], depths[mask]/1e3, 'b', label="t = %.2f Myr" % (_time / 1e6)) # plot
        ax.semilogx(spcrust_segmentations, depths/1e3, 'b', label="t = %.2f Myr" % (_time / 1e6)) # plot
        ax.invert_yaxis()
        ax.set_xlabel("Crust Material in Segments (km^2)")
        ax.set_ylabel("Depth (km)")
        ax.set_xlim([spcrust_segmentations[0]/1e3, spcrust_segmentations[0]*5.0])

    
    def PlotMaterialRate(self, _time, ax, **kwargs):
        '''
        plot rate of tranform
        Inputs:
            ax: an axis is passed in for plotting
        '''
        # need results from two adjacent steps
        dt = kwargs.get('dt', 0.5e6)
        # read data from this step
        data_list, step = self.ha_reader.ExportDataByTime(_time, ["depth", "spcrust", "spharz"])
        depths = data_list[:, 0]
        # t_last: last step to compute the rate of material transformation
        if _time < dt:
            t_last = 0.0
        else:
            t_last = _time - dt
        spcrust_integretions, spcrust_segmentations = self.ha_reader.GetIntegrateArray(_time, "spcrust" ,\
            2, self.geometry_input, self.geometry_extent, Ro=self.Ro)
        spcrust_integretions_last, spcrust_segmentations_last = self.ha_reader.GetIntegrateArray(t_last, "spcrust",\
            2, self.geometry_input, self.geometry_extent, Ro=self.Ro)
        # derive the rate of material transformation
        spcrust_below_depths = spcrust_integretions[-1] - spcrust_integretions
        spcrust_below_depths_last = spcrust_integretions_last[-1] - spcrust_integretions_last
        spcrust_transform_rate = (spcrust_below_depths - spcrust_below_depths_last)/dt
        # plot
        ax.semilogx(spcrust_transform_rate, depths/1e3, 'b', label="t = %.2f Myr" % (_time / 1e6)) # plot
        ax.invert_yaxis()
        ax.set_xlabel("Crust Material transformed (km^2/yr)")
        ax.set_ylabel("Depth (km)")
        ax.set_xlim([spcrust_transform_rate[0]/1e3, spcrust_transform_rate[0]*5.0])
        return step


class PC_MORPH_OPT(PC_OPT_BASE):
    '''
    Define a class to work with json files
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        PC_OPT_BASE.__init__(self)
        self.start = self.number_of_keys()
        self.add_key("time range", list, ["time range"], [], nick="time_range")
        self.add_key("trench position range", list, ["trench position range"], [], nick="tp_range")
        self.add_key("slab depth range", list, ["slab depth range"], [], nick="sd_range")

    def to_init(self):
        '''
        interfaces to the __init__ function
        '''
        case_absolute_paths = self.get_case_absolute_paths()
        return case_absolute_paths

    def to_call(self):
        '''
        interfaces to the __call__ function
        '''
        width = self.values[2]
        output_dir = self.get_output_dir()
        time_range = self.values[self.start]
        tp_range = self.values[self.start + 1]
        sd_range = self.values[self.start + 2]
        return width, output_dir, time_range, tp_range, sd_range


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
        kwargs:
            color_method: use a list of color or the generated values
        '''
        multiple_size = kwargs.get("multiple_size", 1) # get the multiple size factor
        _name = "combine_morphology"
        _title = "Comparing slab morphology results"
        color_method = kwargs.get('color_method', 'list')
        dump_color_to_json = kwargs.get('dump_color_to_json', None)
        save_pdf = kwargs.get('save_pdf', False)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        # initiate
        ni = 3  # number of plots along 1st and 2nd dimension
        nj = 2
        fig, gs, colors = self.initiate_combined_plotting((ni, nj), color_method, dump_color_to_json, multiple_size=multiple_size)
        case_names = []  # names of cases
        for i in range(self.n_cases):
            case_name = os.path.basename(self.cases[i])
            case_names.append(case_name)
        # plot trench position
        ax = fig.add_subplot(gs[1, 0])
        lns = None
        labs = None
        for i in range(self.n_cases):
            if i == 0:
                label = ['Trench Position', 'Slab Depth']
            else:
                label = [None, None]
            case_dir = self.cases[i]
            case_name = os.path.basename(case_dir)
            # plot results and combine
            lns_temp, labs_temp = self.MorphPlotter.PlotTrenchPosition(case_dir, time_range=time_range,\
            tp_range=tp_range, sd_range=sd_range, axis=ax, color=colors[i], label=label)
            if i == 0:
                lns = lns_temp  # record the lables at the start
                labs = labs_temp
            pass
        ax.legend(lns, labs)
        # plot trench velocity
        ax = fig.add_subplot(gs[2, 0])
        lns = None
        labs = None
        for i in range(self.n_cases):
            if i == 0:
                label_all = True
            else:
                label_all = False
            case_dir = self.cases[i]
            case_name = os.path.basename(case_dir)
            # plot results and combine
            self.MorphPlotter.PlotTrenchVelocity(case_dir, time_range=time_range,\
            tp_range=tp_range, sd_range=sd_range, axis=ax, color=colors[i], label_all=label_all)
        ax.legend()
        # plot trench velocity, zoom in
        ax = fig.add_subplot(gs[2, 1])
        for i in range(self.n_cases):
            case_dir = self.cases[i]
            case_name = os.path.basename(case_dir)
            # plot results and combine
            self.MorphPlotter.PlotTrenchVelocity(case_dir, time_range=time_range,\
            tp_range=tp_range, sd_range=sd_range, axis=ax, color=colors[i], label_all=False, fix_v_range=True)
        # plot the color labels
        ax = fig.add_subplot(gs[0, 0])
        PlotColorLabels(ax, case_names, colors)
        # generate figures
        fig_path = os.path.join(output_dir, '%s.png' % _name)
        print("%s: save figure: %s" % (Utilities.func_name(), fig_path))
        plt.savefig(fig_path)
        if save_pdf == True:
            pdf_path = os.path.join(output_dir, '%s.pdf' % _name)
            print("%s: save figure: %s" % (Utilities.func_name(), pdf_path))
            plt.savefig(pdf_path)
        return fig_path

def PlotSlabMaterialTime(case_dir, _time):
    '''
    Plot slab material
    '''
    # mkdir directory
    img_dir = os.path.join(case_dir, 'img')
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    img_sm_dir = os.path.join(img_dir,'slab_material')
    if not os.path.isdir(img_sm_dir):
        os.mkdir(img_sm_dir)
    # initiate plotter
    plotter_material = SLABMATERIAL('slab')
    plotter_material.ReadFile(case_dir)
    fig = plt.figure(tight_layout=True, figsize=(5, 10))
    gs = gridspec.GridSpec(2, 1)
    # plot composition
    ax = fig.add_subplot(gs[0, :])
    step = plotter_material.PlotSlabMaterial(_time, ax)
    ax.legend()
    # plot rate of transformation
    ax = fig.add_subplot(gs[1, :])
    step = plotter_material.PlotMaterialRate(_time, ax)
    ax.legend()
    fileout = os.path.join(img_sm_dir, "s%06d_t%.4e.png" % (step, _time))
    fig.savefig(fileout)
    print("%s: %s generated" % (Utilities.func_name(), fileout))
    pass

def PlotTrenchAgeFromT(case_dir, **kwargs):
    '''
    plot the age of the trench, from the result of the thermal interpretation of the trench temperature
    Inputs:
        kwargs:
            time_interval - interval between steps
    '''
    use_thermal = True  # twik options between thermally interpretation and motion reconstruction
    # let the user check these options
    if use_thermal:
        # use thermal option would fit the temperature at the trench for the individual steps
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        _continue = input("This option will plot the data in the vtk_outputs/trench_T.txt file, \
but will not generarte that file. Make sure all these files are updated, proceed (y/n)?")
        if _continue != 'y':
            print('abort')
            exit(0)
    else:
        _continue = input("This option requires the data in the vtk_outputs/slab_morph.txt file, \
but will not generarte that file. Make sure all these files are updated, proceed (y/n)?")
        if _continue != 'y':
            print('abort')
            exit(0) 
    time_interval = kwargs.get('time_interval', 0.5e6)
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    fig = plt.figure(tight_layout=True, figsize=(5, 5))
    gs = gridspec.GridSpec(1, 1)
    # plot composition
    plotter = SLABPLOT('trench_age_from_T')
    # 0. age
    ax = fig.add_subplot(gs[0, 0])
    plotter.PlotTrenchAge(case_dir, axis=ax, time_interval=time_interval, use_thermal=use_thermal)
    fig_path = os.path.join(img_dir, "trench_age_from_T.png") 
    fig.savefig(fig_path)
    print("%s: save figure %s" % (Utilities.func_name(), fig_path))

def PlotTrenchThermalState(case_dir, **kwargs):
    '''
    plot the age of the trench
    Inputs:
        kwargs:
            time_interval - interval between steps
            silent - function would not ask user input to progress
    '''
    use_thermal = False  # twik options between thermally interpretation and motion reconstruction
    silent = kwargs.get('silent', False)
    # let the user check these options
    if not silent:
        if use_thermal:
            # use thermal option would fit the temperature at the trench for the individual steps
            # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
            _continue = input("This option will plot the data in the vtk_outputs/trench_T.txt file, \
    but will not generarte that file. Make sure all these files are updated, proceed (y/n)?")
            if _continue != 'y':
                print('abort')
                exit(0)
        else:
            _continue = input("This option requires the data in the vtk_outputs/slab_morph.txt file, \
    but will not generarte that file. Make sure all these files are updated, proceed (y/n)?")
            if _continue != 'y':
                print('abort')
                exit(0) 
    time_interval = kwargs.get('time_interval', 0.5e6)
    img_dir = os.path.join(case_dir, "img")
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    fig = plt.figure(tight_layout=True, figsize=(10, 15))
    gs = gridspec.GridSpec(3, 2)
    # plot composition
    plotter = SLABPLOT('trench_thermal_state')
    # 0. age
    ax = fig.add_subplot(gs[0, 0])
    plotter.PlotTrenchAge(case_dir, axis=ax, time_interval=time_interval, use_thermal=use_thermal)
    ax.legend()
    # 1. thermal parameter
    ax = fig.add_subplot(gs[1, 0])
    plotter.PlotThermalParameter(case_dir, axis=ax, time_interval=time_interval, use_thermal=use_thermal)
    ax.legend()
    # 2. thermal parameter, focusing on the stable subduction regem
    ax = fig.add_subplot(gs[2, 0])
    plotter.PlotThermalParameter(case_dir, axis=ax, time_interval=time_interval, use_thermal=use_thermal, time_stable=10e6)
    ax.legend()
    # 3. thermal parameter with the subducting plate velocity
    ax = fig.add_subplot(gs[0, 1])
    plotter.PlotThermalParameter(case_dir, axis=ax, time_interval=time_interval,\
    use_thermal=use_thermal, time_stable=10e6, plot_velocity=True)
    ax.legend()
    # 4. thermal parameter with the dip angle
    ax = fig.add_subplot(gs[1, 1])
    plotter.PlotThermalParameter(case_dir, axis=ax, time_interval=time_interval,\
    use_thermal=use_thermal, time_stable=10e6, plot_dip=True)
    ax.legend()

    fig_path = os.path.join(img_dir, "trench_thermal_state.png") 
    fig.savefig(fig_path)
    print("%s: save figure %s" % (Utilities.func_name(), fig_path))


def PlotMorphAnimeCombined(case_dir, **kwargs):
    '''
    plot slab morphology for making animation
    Inputs:
        case_dir: case directory
        kwargs:
            time_interval - time_interval of plotting
    '''
    # initiate
    time_interval = kwargs.get("time_interval", 5e5)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # SLABPLOT object 
    SlabPlot = SLABPLOT('slab')
    
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    times= Visit_Options.get_times_for_slab_morphology(time_interval=time_interval)
    # print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # debug 
    for _time in times:
        SlabPlot.PlotMorphAnime(case_dir, time=_time)


def PlotTrenchDifferences2dInter1Ma(SlabPlot, case_dir, **kwargs):
    '''
    plot the differences in the trench location since model started (trench migration)
    overlay the curve on an existing axis.
    This function is created for combining results with those from the 3d cases
    '''
    # initiate plot
    _color = kwargs.get('color', "c")
    ax = kwargs.get('axis', None) # for trench position
    ax_twinx = kwargs.get("axis_twinx", None) # for slab depth
    if ax == None:
        raise ValueError("Not implemented")
    # path
    img_dir = os.path.join(case_dir, 'img')
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    morph_dir = os.path.join(img_dir, 'morphology')
    if not os.path.isdir(morph_dir):
        os.mkdir(morph_dir)
    # read inputs
    prm_file = os.path.join(case_dir, 'output', 'original.prm')
    assert(os.access(prm_file, os.R_OK))
    SlabPlot.ReadPrm(prm_file)
    # read parameters
    geometry = SlabPlot.prm['Geometry model']['Model name']
    if geometry == 'chunk':
        Ro = float(SlabPlot.prm['Geometry model']['Chunk']['Chunk outer radius'])
    else:
        Ro = None
    # read data
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph_t1.00e+06.txt')
    assert(os.path.isfile(slab_morph_file))
    SlabPlot.ReadHeader(slab_morph_file)
    SlabPlot.ReadData(slab_morph_file)
    if not SlabPlot.HasData():
        print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
        return 1
    col_pvtu_step = SlabPlot.header['pvtu_step']['col']
    col_pvtu_time = SlabPlot.header['time']['col']
    col_pvtu_trench = SlabPlot.header['trench']['col']
    col_pvtu_slab_depth = SlabPlot.header['slab_depth']['col']
    col_pvtu_sp_v = SlabPlot.header['subducting_plate_velocity']['col']
    col_pvtu_ov_v = SlabPlot.header['overiding_plate_velocity']['col']
    pvtu_steps = SlabPlot.data[:, col_pvtu_step]
    times = SlabPlot.data[:, col_pvtu_time]
    trenches = SlabPlot.data[:, col_pvtu_trench]
    slab_depths = SlabPlot.data[:, col_pvtu_slab_depth]
    time_interval = times[1] - times[0]
    if time_interval < 0.5e6:
        warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)
    if geometry == "chunk":
        trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
    elif geometry == 'box':
        trenches_migration_length = trenches - trenches[0]
    else:
        raise ValueError('Invalid geometry')
    # get_slab_dimensions_2(x, y, Ro, is_chunk)
    ax.plot(times/1e6, trenches_migration_length/1e3, color=_color, label = "2d")
    if ax_twinx is not None:
        ax_twinx.plot(times/1e6, slab_depths/1e3, '--', color=_color)


def PlotSlabDip100km2dInter1Ma(SlabPlot, case_dir, **kwargs):
    '''
    plot the differences in the trench location since model started (trench migration)
    overlay the curve on an existing axis.
    This function is created for combining results with those from the 3d cases
    '''
    # initiate plot
    _color = kwargs.get('color', "c")
    ax = kwargs.get('axis', None) # for trench position
    ax_twinx = kwargs.get("axis_twinx", None) # for slab depth
    if ax == None:
        raise ValueError("Not implemented")
    # path
    img_dir = os.path.join(case_dir, 'img')
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)
    morph_dir = os.path.join(img_dir, 'morphology')
    if not os.path.isdir(morph_dir):
        os.mkdir(morph_dir)
    # read inputs
    prm_file = os.path.join(case_dir, 'output', 'original.prm')
    assert(os.access(prm_file, os.R_OK))
    SlabPlot.ReadPrm(prm_file)
    # read parameters
    geometry = SlabPlot.prm['Geometry model']['Model name']
    # if geometry == 'chunk':
    #     Ro = float(SlabPlot.prm['Geometry model']['Chunk']['Chunk outer radius'])
    # else:
    #     Ro = None
    # read data
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph_t1.00e+06.txt')
    assert(os.path.isfile(slab_morph_file))
    SlabPlot.ReadHeader(slab_morph_file)
    SlabPlot.ReadData(slab_morph_file)
    if not SlabPlot.HasData():
        print("PlotMorph: file %s doesn't contain data" % slab_morph_file)
        return 1
    col_pvtu_step = SlabPlot.header['pvtu_step']['col']
    col_pvtu_time = SlabPlot.header['time']['col']
    col_pvtu_trench = SlabPlot.header['trench']['col']
    col_pvtu_slab_depth = SlabPlot.header['slab_depth']['col']
    # col_pvtu_sp_v = SlabPlot.header['subducting_plate_velocity']['col']
    # col_pvtu_ov_v = SlabPlot.header['overiding_plate_velocity']['col']
    col_100km_dip = SlabPlot.header['100km_dip']['col']
    # pvtu_steps = SlabPlot.data[:, col_pvtu_step]
    times = SlabPlot.data[:, col_pvtu_time]
    trenches = SlabPlot.data[:, col_pvtu_trench]
    # slab_depths = SlabPlot.data[:, col_pvtu_slab_depth]
    dip_100kms = SlabPlot.data[:, col_100km_dip]
    time_interval = times[1] - times[0]
    if time_interval < 0.5e6:
        warnings.warn("Time intervals smaller than 0.5e6 may cause vabriation in the velocity (get %.4e)" % time_interval)

    # Apply a univariate spline to smooth the dip angles
    spline = UnivariateSpline(times / 1e6, dip_100kms, s=0)  # s=0 means interpolation without smoothing
    times_new = np.linspace(times.min(), times.max(), 1000)
    dips_splined = spline(times_new / 1e6)
    # get_slab_dimensions_2(x, y, Ro, is_chunk)
    ax.plot(times/1e6, dip_100kms * 180.0 / np.pi, color=_color, label = "2d")
    ax.plot(times_new/1e6, dips_splined * 180.0 / np.pi, "-.", color=_color, label = "2d")
    # if ax_twinx is not None:
    #     ax_twinx.plot(times/1e6, slab_depths/1e3, '--', color=_color)


def GetSlabDipAt660(case_dir, **kwargs):
    '''
    Get the slab dip angle when reaching 660
    '''
    IndexByValue = lambda array_1d, val: np.argmin(abs(array_1d - val))
    Resample1d = lambda array_1d, n: array_1d[np.ix_(range(0, array_1d.size, n))]

    query_depth = kwargs.get("query_depth", 660e3)
    dip_angle_depth_lookup = kwargs.get("dip_angle_depth_lookup", 660e3)
    dip_angle_depth_lookup_interval = kwargs.get("dip_angle_depth_lookup_interval", 60e3)
    
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret() 
    
    slab_morph_path = os.path.join(case_dir, "vtk_outputs", "slab_morph_t1.00e+05.txt")
    Utilities.my_assert(os.path.isfile(slab_morph_path), SLABPLOT.SlabMorphFileNotExistError, "File %s doesn't exist" % slab_morph_path)
    
    data = np.loadtxt(slab_morph_path)
    steps = data[:, 1]
    times = data[:, 2]
    trenches = data[:, 3]
    slab_depths = data[:, 4]
    
    # time of slab tip reaching _depth km and the index in the list
    sfunc = interp1d(slab_depths, times, assume_sorted=True)
    t_depth = sfunc(query_depth)
    i_depth = IndexByValue(times, t_depth)
    step_depth = steps[i_depth]

    # figure out the snapshot to analyze 
    available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=0.1e6)
    available_pvtu_times = Visit_Options.get_times_for_slab_morphology(time_interval=0.1e6)
    # available_pvtu_times, available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology_outputs(time_interval=0.1e6)
    id = IndexByValue(available_pvtu_times, t_depth)
    vtu_snapshot = available_pvtu_snapshots[id]

    # get the dip angle at _depth km 
    vtu_step, outputs = SlabMorphology_dual_mdd(case_dir, vtu_snapshot, dip_angle_depth_lookup=dip_angle_depth_lookup, dip_angle_depth_lookup_interval=dip_angle_depth_lookup_interval)
    o_list = []
    for entry in outputs.split(' '):
        if entry not in ["", "\n"]:
            o_list.append(entry)
    dip_depth = float(o_list[-1])
    return dip_depth


def get_slab_dimensions_2(x, y, Ro, is_chunk):
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
        r, th1, ph1 = Utilities.cart2sph(x, y, 0.0)
        w = 0.0
        l = Ro * ph1  # Calculate length using the spherical angle ph1
    else:
        # For non-chunk geometry, use Cartesian coordinates directly
        r = y
        w = 0.0
        l = x 

    return r, w, l


# todo_o_env
def SlabEnvelopRetrivePoints(local_dir: str, _time: float, Visit_Options: object, depths: float, **kwargs: dict) -> tuple:
    '''
    Retrieves the point of the slab envelop at give depth

    Parameters:
        local_dir (str): Directory where output files are located.
        _time (float): The time at which to retrieve the heat flow profile.
        Visit_Options (object): An object containing options for retrieving data.
        depth: the depth of the point to retrive
        kwargs (dict): Additional keyword arguments; accepts 'phi_diff' (float, default: 5.0).

    Returns:
        tuple: A tuple containing two masked arrays: heat fluxes and corresponding phi values.
    '''
    Ro = Visit_Options.options["OUTER_RADIUS"]

    _time1, timestep, vtu_step = Visit_Options.get_timestep_by_time(_time)
    filein = os.path.join(local_dir, "vtk_outputs", "slab_env_%05d.txt" % vtu_step)
    Utilities.my_assert(os.path.isfile(filein), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), filein))

    # Retrieve slab envelops 
    # and interpolate the points based on the give depths
    data = np.loadtxt(filein)
    X1 = data[:, 2]
    Y1 = data[:, 3]

    # L1: box geometry - X dimension; chunk geometry - phi dimension
    if Visit_Options.options["GEOMETRY"] == "box":
        Ys = np.interp(Ro-depths, Y1,X1) 
        return Ys
    elif Visit_Options.options["GEOMETRY"] == "chunk":
        R0, _, Phi0 = Utilities.cart2sph(X1,Y1,np.zeros(X1.shape))
        Phis = np.interp(depths, Ro-R0, Phi0)
        return Phis
    else:
        raise NotImplementedError

# todo_shallow 
def minimum_distance_array(a, x0, y0, z0):
    """
    Calculate the minimum distance from a reference point (x0, y0, z0) 
    to the points in array 'a', where 'a' is a numpy array of shape (n, 3).
    
    Parameters:
        a (numpy.ndarray): Array of shape (n, 3), where each row is [x, y, z] coordinates.
        x0 (float): x-coordinate of the reference point.
        y0 (float): y-coordinate of the reference point.
        z0 (float): z-coordinate of the reference point.
    
    Returns:
        float: The minimum distance from (x0, y0, z0) to the points in 'a'.
    """
    # Calculate the squared distances to avoid unnecessary square roots
    squared_distances = (a[:, 0] - x0)**2 + (a[:, 1] - y0)**2 + (a[:, 2] - z0)**2
    
    # Find the minimum squared distance and take its square root to get the actual distance
    min_distance = np.sqrt(np.min(squared_distances))
    min_index = np.argmin(squared_distances)

    return min_index, min_distance


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
    parser.add_argument('-vs', '--vtu_step', type=int,
                        default=0,
                        help='vtu_step')
    parser.add_argument('-vss', '--vtu_snapshot', type=int,
                        default=0,
                        help='vtu_snapshot')
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='A json file for configuration')
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=1e6,
                        help='Time interval, affecting the time steps to visualize')
    parser.add_argument('-te', '--time_end', type=float,
                        default=60e6,
                        help='Time end, affecting the time steps to visualize')
    parser.add_argument('-ts', '--time_start', type=float,
                        default=0.0,
                        help='Time end, affecting the time steps to visualize')
    parser.add_argument('-t', '--time', type=float,
                        default=0.0,
                        help='Time')
    parser.add_argument('-sp', '--save_pdf', type=int,
                        default=0,
                        help='save file to pdf format')
    parser.add_argument('-fmdt', '--findmdd_tolerance', type=float,
                        default=0.05,
                        help='Tolerence when finding the mdd')
    parser.add_argument('-mdx0', '--mdd_dx0', type=float,
                        default=10e3,
                        help='Distance from surface when finding the mdd (smaller x side)')
    parser.add_argument('-mdx1', '--mdd_dx1', type=float,
                        default=10e3,
                        help='Distance from surface when finding the mdd (larger x side)')
    # todo_T
    parser.add_argument('-rp', '--run_parallel', type=int,
                        default=1,
                        help='run in parallel')
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
        # example:
        Visit_Options = VISIT_OPTIONS(arg.inputs)
        # call function
        Visit_Options.Interpret()
        ofile = arg.outputs
        SlabAnalysis(arg.inputs, arg.vtu_snapshot, ofile)
    elif _commend == 'plot_slab_envelops': 
        # plot slab envelops
        PlotSlabEnvelops(arg.inputs, arg.outputs, arg.vtu_step, include_internal=True)
    elif _commend == 'plot_slab_forces':
        # plot slab forces
        PlotSlabForces(arg.inputs, arg.outputs)
    elif _commend == "plot_slab_case_step":
        PlotSlabForcesCase(arg.inputs, arg.vtu_step, output_slab=True)
    elif _commend == "plot_slab_shape":
        PlotSlabShape(arg.inputs, arg.vtu_step)
    elif _commend == 'morph_step':
        # slab_morphology, input is the case name
        SlabMorphology_dual_mdd(arg.inputs, int(arg.vtu_snapshot), rewrite=1, findmdd=True, project_velocity=True,\
            findmdd_tolerance=arg.findmdd_tolerance, mdd_dx0=arg.mdd_dx0, mdd_dx1=arg.mdd_dx1, output_ov_ath_profile=True, output_slab='txt')
    # todo_env
    elif _commend == 'morph_case':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=1, findmdd=True, time_interval=arg.time_interval, project_velocity=True,\
            findmdd_tolerance=arg.findmdd_tolerance, output_ov_ath_profile=True, output_slab='txt', find_shallow_trench=True)
    elif _commend == 'morph_case_parallel':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=1, findmdd=True, time_interval=arg.time_interval, project_velocity=True,\
            use_parallel=True, file_tag='interval', findmdd_tolerance=arg.findmdd_tolerance, output_ov_ath_profile=True, output_slab='txt', find_shallow_trench=True)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotMorph(arg.inputs, save_pdf=arg.save_pdf)
    elif _commend == 'plot_morph_animation':
        PlotMorphAnimeCombined(arg.inputs)
    elif _commend == 'plot_morph_publication':
        # plot slab morphology for animation
        plt.style.use('publication_2d_morph')
        SlabPlot = SLABPLOT('slab')
        time_interval = 0.5e6
        time_range = [0.0, 44.6e6]
        time_markers = [2.7e6, 13.0e6, 20.3e6]
        vlim = [-10.0, 20.0]
        SlabPlot.PlotMorphPublication(arg.inputs, save_pdf=True, time_interval=time_interval, time_range=time_range, time_markers=time_markers, vlim=vlim)
    elif _commend == 'plot_wedge_T':
        PlotWedgeTCase(arg.inputs, time_interval=arg.time_interval)
    elif _commend == 'combine_slab_morph':
        # combine plot of slab morphology
        _continue = input("This option will plot the data in the vtk_outputs/slab_morph.txt file, \
but will not generarte that file. Make sure all these files are updated, proceed (y/n)?")
        if _continue == 'y':
            PlotCombineExecute(PLOT_COMBINE_SLAB_MORPH, PC_MORPH_OPT, "slab_morph", Utilities.re_neat_word(arg.json))
        else:
            print("abort")
    elif _commend == "plot_shear_zone":
        ShearZoneGeometry(arg.inputs, int(arg.vtu_snapshot))
    elif _commend == "plot_shear_zone_case":
        ShearZoneGeometryCase(arg.inputs, time_interval=arg.time_interval, time_end=arg.time_end, time_start=arg.time_start)
    elif _commend == "plot_slab_material_time":
        PlotSlabMaterialTime(arg.inputs, arg.time)
    elif _commend == "plot_slab_temperature":
        PlotSlabTemperature(arg.inputs, int(arg.vtu_snapshot))
    elif _commend == "slab_temperature_case":
        # todo_T
        SlabTemperatureCase(arg.inputs, rewrite=1, time_interval=arg.time_interval, offsets=[-5e3, -10e3], run_parallel=arg.run_parallel)
    elif _commend == "mantle_wedge_T_case":
        WedgeTCase(arg.inputs, time_interval=arg.time_interval)
    elif _commend == "plot_slab_temperature_case":
        if arg.time_start > 0 or arg.time_end < 60e6:
            PlotSlabTemperatureCase(arg.inputs, time_range=[arg.time_start, arg.time_end], plot_eclogite=True)
        else:
            PlotSlabTemperatureCase(arg.inputs)
    elif _commend == "trench_T_case":
        TrenchTCase(arg.inputs, time_interval=arg.time_interval)
    elif _commend == "plot_slab_age_from_T":
        PlotTrenchAgeFromT(arg.inputs, time_interval=arg.time_interval)
    elif _commend == "plot_slab_thermal":
        PlotTrenchThermalState(arg.inputs, time_interval=arg.time_interval)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
