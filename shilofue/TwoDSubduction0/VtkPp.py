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
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import vtk
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import colors as mcolors
import shilofue.VtkPp as VtkPp
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from numpy import linalg as LA 
import multiprocessing
#### self
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile, ParseFromDealiiInput
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
  - python -m shilofue.TwoDSubduction0.VtkPp morph_step -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 -vss 105 \n\
\n\
  - Plot the morphology of the slab: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_morph -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 \n\
\n\
        ")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
        slab_cells: cell id of internal points in the slab
        slab_envelop_cell_list0: cell id of slab envelop with smaller theta, slab bottom
        slab_envelop_cell_list1: cell id of slab envelop with bigger theta, slab surface
        slab_trench: trench position, theta for a 'chunk' model and x for a 'box' model
        coord_100: position where slab is 100km deep, theta for a 'chunk' model and x for a 'box' model
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
        self.surface_cells = []
        self.slab_envelop_cell_list0 = []
        self.slab_envelop_cell_list1 = []
        self.slab_depth = None
        self.slab_trench = None
        self.coord_100 = None 
        self.dip_100 = None
        self.slab_shallow_cutoff = 50e3  # depth limit to slab
        self.slab_envelop_interval = 5e3
        default_gravity_file = os.path.join(Utilities.var_subs('${ASPECT_SOURCE_DIR}'),\
        "data", "gravity-model", "prem.txt") 
        gravity_file = kwargs.get('gravity_file', default_gravity_file)
        assert(os.path.isfile(gravity_file))
        self.ImportGravityData(gravity_file)

    def PrepareSlab(self, slab_field_names, **kwargs):
        '''
        prepare slab composition
        '''
        assert(self.include_cell_center)
        slab_threshold = kwargs.get('slab_threshold', 0.2)
        points = vtk_to_numpy(self.i_poly_data.GetPoints().GetData())
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())
        point_data = self.i_poly_data.GetPointData()
        cell_point_data = self.c_poly_data.GetPointData()
        # slab composition field
        slab_field = VtkPp.OperateDataArrays(cell_point_data, slab_field_names,\
        [0 for i in range(len(slab_field_names) - 1)])
        # add cells by composition
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

    
    def ExportSlabInfo(self):
        '''
        Output slab information
        '''
        return self.trench, self.slab_depth, self.dip_100


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
        indent = kwargs.get('indent', 0)
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


    def SlabBuoyancy(self, v_profile, depth_increment):
        '''
        slab buoyancy
        v_profile: vertical profile containing the reference profile
        rs : radius segments for computing buoyancy
        '''
        grav_acc = 10.0
        assert(self.include_cell_center)
        assert(len(self.slab_cells) > 0)
        n_depth = int(np.ceil(self.slab_depth / depth_increment))
        buoyancies = np.zeros(n_depth)
        depths = []  # construct depth array
        for i in range(n_depth):
            depth = (i + 0.5) * depth_increment
            depths.append(depth)
        depths = np.array(depths)
        centers = vtk_to_numpy(self.c_poly_data.GetPoints().GetData())  # note these are data mapped to cell center
        density_data = vtk_to_numpy(self.c_poly_data.GetPointData().GetArray('density'))
        density_ref_func = v_profile.GetFunction('density')
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

####
# Utilities functions
####
def get_r(x, y, geometry): 
    '''
    Get r (the first coordinate)
    Inputs:
        x - x coordinate
        y - y coordinate
        geometry - 'chunk' or 'box'
    '''
    if geometry == 'chunk':
        r = (x*x + y*y)**0.5
    elif geometry == 'box':
        r = y
    else:
        raise ValueError("not implemented")
    return r


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
    '''
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    assert(os.path.isfile(filein))
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    # vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    VtkP = VTKP(geometry=geometry, Ro=Ro)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    # generate outputs
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e%-14.4e\n"\
    % (vtu_step, step, _time, trench, slab_depth, dip_100)
    print(outputs) # debug
    return vtu_step, outputs


def SlabAnalysis(case_dir, vtu_snapshot, o_file, **kwargs):
    '''
    Perform analysis on the slab, this would output a file including the
    buoyancy forces of the slab and the pressures on the slab surface.
    Inputs:
        kwargs(dict):
            output_slab - output slab file
    '''
    # assert something
    indent = kwargs.get("indent", 0)  # indentation for outputs
    print("%s%s: Start" % (indent*" ", Utilities.func_name()))
    output_slab = kwargs.get('output_slab', False)
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
    # initiate class
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'p', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # output slab profile
    if output_slab:
        slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
        slab_internal = VtkP.ExportSlabInternal(output_xy=True)
        o_slab_env0 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env0_%05d.txt" % (vtu_step)) # envelop 0
        o_slab_env1 = os.path.join(case_dir,\
            "vtk_outputs", "slab_env1_%05d.txt" % (vtu_step)) # envelop 1
        o_slab_in = os.path.join(case_dir,\
            "vtk_outputs", "slab_internal_%05d.txt" % (vtu_step)) # envelop 1
        np.savetxt(o_slab_env0, slab_envelop0)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env0))
        np.savetxt(o_slab_env1, slab_envelop1)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_env1))
        np.savetxt(o_slab_in, slab_internal)
        print("%s%s: write file %s" % (indent*" ", Utilities.func_name(), o_slab_in))
    # buoyancy
    r0_range = [6371e3 - 2890e3, 6371e3]
    Ro = 6371e3
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    total_buoyancy, b_profile = VtkP.SlabBuoyancy(v_profile, 5e3)  # test 5e3, 50e3
    depths_o = b_profile[:, 0]  # use these depths to generate outputs
    buoyancies = b_profile[:, 1]
    buoyancy_gradients = buoyancies / (depths_o[1] - depths_o[0])  # gradient of buoyancy
    # pressure 
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()  # raw data on the envelop and output
    fileout = os.path.join(output_path, 'slab_pressures0_%05d.txt' % (vtu_step))
    depths0, thetas0, ps0= SlabPressures(VtkP, slab_envelop0, fileout=fileout, indent=4)  # depth, dip angle and pressure
    fileout = os.path.join(output_path, 'slab_pressures1_%05d.txt' % (vtu_step))
    depths1, thetas1, ps1 = SlabPressures(VtkP, slab_envelop1, fileout=fileout, indent=4)
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
    fileout = kwargs.get('fileout', None)  # debug
    indent = kwargs.get('indent', 0)
    rs_n = 5 # resample interval
    ip_interval = 1e3  # interval for interpolation
    # resample the original envelop dataset
    n_point = slab_envelop.shape[0]
    rs_idx = range(0, n_point, rs_n)
    slab_envelop_rs = slab_envelop[np.ix_(rs_idx, [0, 1])]
    slab_env_polydata = VtkPp.InterpolateGrid(VtkP.i_poly_data, slab_envelop_rs, quiet=True) # note here VtkPp is module shilofue/VtkPp, while the VtkP is the class
    temp_vtk_array = slab_env_polydata.GetPointData().GetArray('p')
    env_ps  = vtk_to_numpy(temp_vtk_array)
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
        p_d = p - p_static # dynamic pressure
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
    # get all available snapshots
    # the interval is choosen so there is no high frequency noises
    time_interval_for_slab_morphology = 0.5e6
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    available_pvtu_steps = [i - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) for i in available_pvtu_snapshots]
    # get where previous session ends
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
    # Initiation Wrapper class for parallel computation
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_morph', SlabMorphology, if_rewrite=True)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file
    if if_rewrite:
        if os.path.isfile(slab_morph_file):
            print("%s: Delete old slab_morph.txt file." % Utilities.func_name())
            os.remove(slab_morph_file)  # delete slab morph file
        ParallelWrapper.delete_temp_files(available_pvtu_steps)  # delete intermediate file if rewrite
    num_cores = multiprocessing.cpu_count()
    # loop for all the steps to plot
    # Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_step)\
    # for pvtu_step in available_pvtu_steps)  # first run in parallel and get stepwise output
    ParallelWrapper.clear()
    for pvtu_step in available_pvtu_steps:  # then run in on cpu to assemble these results
        ParallelWrapper(pvtu_step)
    pvtu_steps_o, outputs = ParallelWrapper.assemble()
    # last, output
    # header
    file_header = "# 1: pvtu_step\n# 2: step\n# 3: time (yr)\n# 4: trench (rad)\n# 5: slab depth (m)\n# 6: 100km dip (rad)\n"
    output_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
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
        self.wedge_T_reader = LINEARPLOT('wedge_T')

    
    def ReadWedgeT(self, case_dir, min_pvtu_step, max_pvtu_step, **kwargs):
        time_interval_for_slab_morphology = 0.5e6  # hard in
        i = 0
        initial_adaptive_refinement = int(self.prm['Mesh refinement']['Initial adaptive refinement'])
        geometry = self.prm['Geometry model']['Model name']
        if geometry == 'chunk':
            Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        elif geometry == 'box':
            Do = float(self.prm['Geometry model']['Box']['Y extent'])
        else:
            raise ValueError('Invalid geometry')
        Visit_Options = VISIT_OPTIONS(case_dir)
        Visit_Options.Interpret()
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        # for pvtu_step in range(min_pvtu_step + initial_adaptive_refinement, max_pvtu_step + initial_adaptive_refinement + 1):
        for pvtu_step in available_pvtu_snapshots:
            file_in_path = os.path.join(case_dir, 'vtk_outputs', 'wedge_T100_%05d.txt' % pvtu_step)
            Utilities.my_assert(os.access(file_in_path, os.R_OK), FileExistsError, "File %s doesn\'t exist" % file_in_path)
            self.wedge_T_reader.ReadHeader(file_in_path)
            self.wedge_T_reader.ReadData(file_in_path)
            col_x = self.wedge_T_reader.header['x']['col']
            col_y = self.wedge_T_reader.header['y']['col']
            col_T = self.wedge_T_reader.header['T']['col']
            xs = self.wedge_T_reader.data[:, col_x]
            ys = self.wedge_T_reader.data[:, col_y]
            if i == 0: 
                rs = (xs**2.0 + ys**2.0)**0.5
                if geometry == 'chunk':
                    depthes = Ro - rs # compute depth
                elif geometry == 'box':
                    depthes = Do - rs # compute depth
                else:
                    raise ValueError('Invalid geometry')
                # Ts = np.zeros((depthes.size, max_pvtu_step - min_pvtu_step + 1))
                Ts = np.zeros((depthes.size, len(available_pvtu_snapshots)))
            Ts[:, i] = self.wedge_T_reader.data[:, col_T]
            i += 1
        return depthes, Ts

    def PlotMorph(self, case_dir, **kwargs):
        '''
        Inputs:
            case_dir (str): directory of case
        kwargs(dict):
            defined but not used
        '''
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
        # trench velocity
        # start figure
        fig = plt.figure(tight_layout=True, figsize=(15, 10)) 
        fig.subplots_adjust(hspace=0)
        gs = gridspec.GridSpec(3, 2) 
        # 1: trench & slab movement
        ax = fig.add_subplot(gs[0, 0:2]) 
        ax_tx = ax.twinx()
        _color = 'tab:blue'
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color=_color, label='trench position (km)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylabel('Trench Position (km)', color=_color)
        ax.tick_params(axis='x', labelbottom=False) # labels along the bottom edge are off
        ax.tick_params(axis='y', labelcolor=_color)
        ax.grid()
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, 'k-', label='slab depth (km)')
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        # ax.legend(lns, labs)
        # 2: velocity
        ax = fig.add_subplot(gs[1, 0:2]) 
        _color = 'tab:blue'
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color=_color, label='trench velocity (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='trench velocity (cm/yr)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylim((-10, 10))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.set_xlabel('Times (Myr)')
        ax.grid()
        # 2.1: velocity smaller, no y limit, to show the whole curve
        ax = fig.add_subplot(gs[2, 0]) 
        _color = 'tab:blue'
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color=_color, label='trench velocity (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='trench velocity (cm/yr)')
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_ylabel('Velocity (whole, cm/yr)')
        ax.grid()
        # 3: wedge temperature
#        depthes, Ts = self.ReadWedgeT(case_dir, int(pvtu_steps[0]), int(pvtu_steps[-1]))
#        tt, dd = np.meshgrid(times, depthes)
#        ax = fig.add_subplot(gs[2, 0:2]) 
#        h = ax.pcolormesh(tt/1e6,dd/1e3,Ts, shading='gouraud') 
#        ax.invert_yaxis()
#        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
#        ax.set_xlabel('Times (Myr)')
#        ax.set_ylabel('Depth (km)')
#        ax = fig.add_subplot(gs[2, 2])
#        ax.axis('off') 
#        fig.colorbar(h, ax=ax, label='T (K)') 
        fig.tight_layout()
        # save figure
        o_path = os.path.join(morph_dir, 'trench.png')
        plt.savefig(o_path)
        print("%s: figure %s generated" % (Utilities.func_name(), o_path))


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
        SlabMorphology(arg.inputs, int(arg.vtu_snapshot), rewrite=1)
    elif _commend == 'morph_case':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotMorph(arg.inputs)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()