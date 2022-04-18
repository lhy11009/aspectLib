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
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import numpy as np
import vtk
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.VtkPp as VtkPp
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
from numpy import linalg as LA 

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
            -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0 -s 100 -o foo.txt\n\
\n\
  - plot the slab forces: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_forces -i foo.txt -o foo.png\n\
\n\
  - Perform analysis on a single step for one case: \n\
        python -m shilofue.TwoDSubduction0.VtkPp plot_slab_case_step -i \n\
            /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT3/eba_cdpt_SA80.0_OA40.0  -s 100\n\
\n\
        ")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
        slab_cells: cell id of internal points in the slab
        slab_envelop_cell_list0: cell id of slab envelop with smaller theta
        slab_envelop_cell_list1: cell id of slab envelop with bigger theta
    '''
    def __init__(self):
        VtkPp.VTKP.__init__(self)
        self.slab_cells = []
        self.surface_cells = []
        self.slab_envelop_cell_list0 = []
        self.slab_envelop_cell_list1 = []
        self.slab_depth = None
        self.slab_shallow_cutoff = 50e3  # depth limit to slab
        self.slab_envelop_interval = 5e3
        self.Ro = 6371e3

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
        is_first = True
        for field_name in slab_field_names:
            if is_first:
                slab_field = vtk_to_numpy(cell_point_data.GetArray(field_name))
                is_first = False
            else:
                slab_field += vtk_to_numpy(cell_point_data.GetArray(field_name))
        Ts = self.i_poly_data.GetPointData().GetArray("T")  # temperature field
        # add cells by composition
        min_r = self.Ro
        for i in range(self.i_poly_data.GetNumberOfCells()):
            cell = self.i_poly_data.GetCell(i)
            id_list = cell.GetPointIds()  # list of point ids in this cell
            x = centers[i][0]
            y = centers[i][1]
            r = (x*x + y*y)**0.5
            slab = slab_field[i]
            if slab > slab_threshold and (self.Ro - r) > self.slab_shallow_cutoff:
                self.slab_cells.append(i)
                if r < min_r:
                    min_r = r
        self.slab_depth = self.Ro - min_r
        # get slab envelops
        total_en_interval = int((self.slab_depth - self.slab_shallow_cutoff) // self.slab_envelop_interval + 1)
        slab_en_cell_lists = [ [] for i in range(total_en_interval) ]
        for id in self.slab_cells:
            x = centers[id][0]  # first, separate cells into intervals
            y = centers[id][1]
            r = (x*x + y*y)**0.5
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
                theta = np.arctan2(x, y)
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
    
    def ExportSlabEnvelopCoord(self):
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
            r = (x*x + y*y)**0.5
            i_r = int(np.floor((self.Ro - r) / depth_increment))
            density = density_data[i]
            density_ref = density_ref_func(r)
            cell_size = self.cell_sizes[i]  # temp
            buoyancy = grav_acc * (density - density_ref) * cell_size
            buoyancies[i_r] += buoyancy
            total_buoyancy += buoyancy
        b_profile = np.zeros((n_depth, 2))
        b_profile[:, 0] = depths
        b_profile[:, 1] = buoyancies
        return total_buoyancy, b_profile


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
    depth_interval = data[1, 0] - data[0, 0]
    buoyancies = data[:, 1]
    total_buoyancy = LA.norm(buoyancies, 1)  # total buoyancy
    buoyancie_gradients = buoyancies / depth_interval
    v_zeros = np.zeros(data.shape[0])
    fig, ax = plt.subplots()
    ax.plot(buoyancie_gradients, depths/1e3, 'b', label='Buoyancy (N/m2)')
    ax.plot(v_zeros, depths/1e3, 'k--')
    ax.invert_yaxis()
    ax.set_title('Total buoyancy: %.4e N/m2' % total_buoyancy)
    ax.set_xlabel('Force (N/m2)')
    ax.set_ylabel('Depth (km)')
    plt.savefig(fileout)
    print("PlotSlabForces: plot figure", fileout)


def SlabAnalysis(case_dir, vtu_step, o_file, **kwargs):
    '''
    Perform analysis on the slab
    Inputs:
        kwargs(dict):
            output_slab - output slab file
    '''
    # assert something
    output_slab = kwargs.get('output_slab', False)
    output_path = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution",\
         "solution-%05d.pvtu" % (vtu_step))
    assert(os.path.isfile(filein))
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
        print("\t%s: write file %s" % (Utilities.func_name(), o_slab_env0))
        np.savetxt(o_slab_env1, slab_envelop1)
        print("\t%s: write file %s" % (Utilities.func_name(), o_slab_env1))
        np.savetxt(o_slab_in, slab_internal)  # todo
        print("\t%s: write file %s" % (Utilities.func_name(), o_slab_in))
    # buoyancy
    r0_range = [6371e3 - 2890e3, 6371e3]
    Ro = 6371e3
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    total_buoyancy, b_profile = VtkP.SlabBuoyancy(v_profile, 5e3)
    header = "# depth (m)\n# buoyancy (N/m)\n"
    with open(o_file, 'w') as fout:
        fout.write(header)  # output header
    with open(o_file, 'a') as fout:
        np.savetxt(fout, b_profile, fmt="%20.8e")  # output data
    print("%s: write file %s" % (Utilities.func_name(), o_file))
    # pressure 
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord()
    fileout = os.path.join(output_path, 'slab_pressures0_%05d.txt' % (vtu_step))
    SlabPressures(VtkP, slab_envelop0, fileout=fileout, indent=4)
    fileout = os.path.join(output_path, 'slab_pressures1_%05d.txt' % (vtu_step))
    SlabPressures(VtkP, slab_envelop1, fileout=fileout, indent=4)


def SlabPressures(VtkP, slab_envelop, **kwargs):
    '''
    extract slab pressures
    Inputs:
        VtkP: VTKP class
        slab_envelop: slab envelop coordinates (x and y)
    '''
    # pressure gradient
    Ro = 6371e3
    fileout = kwargs.get('fileout', None)  # debug
    indent = kwargs.get('indent', 0)
    slab_env_polydata = VtkPp.InterpolateGrid(VtkP.i_poly_data, slab_envelop)
    temp_vtk_array = slab_env_polydata.GetPointData().GetArray('p')
    env_ps  = vtk_to_numpy(temp_vtk_array)
    temp_vtk_array = slab_env_polydata.GetPointData().GetArray('p')
    # env1_ps = vtk_to_numpy(temp_vtk_array.GetData())
    depths = np.zeros(slab_envelop.shape[0]) # data on envelop0
    ps = np.zeros((slab_envelop.shape[0], 3)) # pressure, horizontal & vertical components
    thetas = np.zeros((slab_envelop.shape[0], 1))
    is_first = True
    for i in range(0, slab_envelop.shape[0]):
        x = slab_envelop[i, 0]
        y = slab_envelop[i, 1]
        depth = Ro - (x * x + y * y)**0.5  # depth of this point
        p = env_ps[i]  # pressure of this point
        depths[i] = depth
        dx = 0.0  # coordinate differences
        dy = 0.0
        if is_first:
            xnext = slab_envelop[i+1, 0]  # coordinates of this and the last point
            ynext = slab_envelop[i+1, 1]
            dx = xnext - x
            dy = ynext - y
            is_first = False
        else: 
            xlast = slab_envelop[i-1, 0]  # coordinates of this and the last point
            ylast = slab_envelop[i-1, 1]
            dx = x - xlast
            dy = y - ylast
        theta = np.arctan2(dx, dy)
        thetas[i, 0] = theta
        p_v = p * np.cos(theta)
        p_h = p * np.sin(theta)
        ps[i, 0] = p
        ps[i, 1] = p_h
        ps[i, 2] = p_v
    temp = np.concatenate((slab_envelop, thetas), axis=1)
    data_env0 = np.concatenate((temp, ps), axis=1)  # assemble all the data
    c_out = data_env0.shape[1]
    start = np.ceil(depths[0]/1e3) * 1e3
    end = np.floor(depths[-1]/1e3) * 1e3
    n_out = int((end-start) / 1e3)
    data_env0_out = np.zeros((n_out, c_out+1))
    depths_out = np.arange(start, end, 1e3)
    data_env0_out[:, 0] = depths_out
    for j in range(c_out):
        data_env0_out[:, j+1] = np.interp(depths_out, depths, data_env0[:, j]) # interpolate to regular grid
    if fileout != None:
        header = "# 1: depth (m)\n# 2: x (m)\n# 3: y (m)\n# 4: theta_v \n# 5: p (Pa)\n# 6: p_h (Pa) \n# 7: p_v (Pa)\n"
        with open(fileout, 'w') as fout:
            fout.write(header)
        with open(fileout, 'a') as fout: 
            np.savetxt(fout, data_env0_out, fmt='%.4e\t')
        print("%s%s: write output %s" % (' '*indent, Utilities.func_name(), fileout))
    


def PlotSlabCase(case_dir, step, **kwargs):
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
    step = step
    vtu_step = int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) + step
    vtk_output_dir = os.path.join(case_dir, 'vtk_outputs')
    if not os.path.isdir(vtk_output_dir):
        os.mkdir(vtk_output_dir)
    ofile = os.path.join(vtk_output_dir, "slab_forces_%05d" % vtu_step)
    SlabAnalysis(case_dir, vtu_step, ofile, output_slab=output_slab)
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
        step = arg.step
        vtu_step = int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) + step
        ofile = arg.outputs
        SlabAnalysis(arg.inputs, vtu_step, ofile)
    elif _commend == 'plot_slab_forces':
        # plot slab forces
        PlotSlabForces(arg.inputs, arg.outputs)
    elif _commend == "plot_slab_case_step":
        PlotSlabCase(arg.inputs, arg.step, output_slab=True)
    elif _commend == "plot_slab_shape":
        PlotSlabShape(arg.inputs, arg.step)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()