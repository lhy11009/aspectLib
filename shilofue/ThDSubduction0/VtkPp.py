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
from matplotlib import pyplot as plt
from matplotlib import cm, gridspec
from numpy import linalg as LA 
import multiprocessing
import time
#### self
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.ThDSubduction0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile
from shilofue.Plot import LINEARPLOT
from shilofue.PlotCombine import PLOT_COMBINE
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
    plot trench positions, assuming the trench.txt files have already been generated: \n\
        the -ti option gives a time interval\n\
\n\
        python -m shilofue.ThDSubduction0.VtkPp plot_trench -i /mnt/lochy0/ASPECT_DATA/ThDSubduction/gmg_test_stampede/test_ThD_gmg_mv1e20_picard_correction_side_plate -ti 2e6 \n\
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
        if i > 0:
            outputs += "\n"
        outputs += str(trench_coords_x[i]) + " " + str(trench_coords_y[i]) + " " + str(trench_coords_z[i])
    fileout = os.path.join(output_path, "trench_%05d.txt" % vtu_snapshot)
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    print("File output for trench positions: %s" % fileout)
    end = time.time()
    print("Write trench positions, takes %.2f s" % (end - start))
    start = end
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
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_morph', SlabMorphology, if_rewrite=True)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file

#    if if_rewrite:
#        if os.path.isfile(slab_morph_file):
#            print("%s: Delete old slab_morph.txt file." % Utilities.func_name())
#            os.remove(slab_morph_file)  # delete slab morph file
#        ParallelWrapper.delete_temp_files(available_pvtu_snapshots)  # delete intermediate file if rewrite

    num_cores = multiprocessing.cpu_count()
    # loop for all the steps to plot
    # Parallel(n_jobs=num_cores)(delayed(ParallelWrapper)(pvtu_step)\
    # for pvtu_step in available_pvtu_steps)  # first run in parallel and get stepwise output
    ParallelWrapper.clear()
    for pvtu_snapshot in available_pvtu_snapshots:  # then run in on cpu to assemble these results
        ParallelWrapper(pvtu_snapshot)
        ParallelWrapper(pvtu_snapshot + step_for_derivatives)  # for computing the derivatives, might be bugs

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
        # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
        available_pvtu_snapshots = Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
        n_snapshots = len(available_pvtu_snapshots)
        print("available_pvtu_snapshots: ", available_pvtu_snapshots)  # output the availabe pvtu snapshots
        fig = plt.figure(tight_layout=True, figsize=(5, 10)) 
        gs = gridspec.GridSpec(2, 1) 
        # 1. trench position
        ax = fig.add_subplot(gs[0, 0]) 
        normalizer = [ float(i)/(n_snapshots-1) for i in range(n_snapshots)] 
        colors = cm.rainbow(normalizer) 
        i = 0
        for vtu_snapshot in available_pvtu_snapshots:
            # plot each snapsh
            _time = time_interval_for_slab_morphology * i
            _color = colors[i]
            ax = self.PlotTrenchPositionStep(case_dir, vtu_snapshot, axis=ax, color=_color, label="Time = %.2f Myr" % (_time / 1e6))
            i += 1
        ax.legend()
        ax.set_xlim([trench_initial_position/1e3 - 1000, trench_initial_position/1e3 + 1000]) # set the x axis to center around the inital trench
        # 2. trench velocity
        ax = fig.add_subplot(gs[1, 0]) 
        ax = self.PlotTrenchVelocity(case_dir, time_interval_for_slab_morphology, axis=ax)
        # save image
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_img_dir = os.path.join(img_dir, "morphology")
        if not os.path.isdir(morph_img_dir):
            os.mkdir(morph_img_dir)
        fileout = os.path.join(morph_img_dir, "trench.png")
        fig.savefig(fileout)
        print("image output: ", fileout)


    def PlotTrenchPositionStep(self, case_dir, vtu_snapshot, **kwargs):
        '''
        plot trench position for a single step
        '''
        _color = kwargs.get("color", None)
        _label = kwargs.get("label", None)
        filein = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % vtu_snapshot)
        assert(os.path.isfile(filein))
        # initiate
        ax = kwargs.get('axis', None)
        if ax == None:
            raise ValueError("Not implemented")
        data = np.loadtxt(filein)
        xs = data[:, 0]
        ys = data[:, 1]
        ax.plot(xs/1e3, ys/1e3, '.', color=_color, label=_label)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        return ax
    
    
    def PlotTrenchVelocity(self, case_dir, time_interval_for_slab_morphology, **kwargs):
        '''
        plot trench position for a single step
        '''
        Myr = 1e6
        cm = 1e-2
        _color = kwargs.get("color", None)
        _label = kwargs.get("label", None)
        step_for_derivatives = kwargs.get('step_for_derivatives', 5)
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
        for vtu_snapshot in available_pvtu_snapshots: 
            filein0 = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % vtu_snapshot)
            filein1 = os.path.join(case_dir, "vtk_outputs", "trench_%05d.txt" % (vtu_snapshot+step_for_derivatives))
            assert(os.path.isfile(filein0))
            assert(os.path.isfile(filein1))
            _time, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot)
            _time1, step = Visit_Options.get_time_and_step_by_snapshot(vtu_snapshot+1)
            dt = _time1 - _time
            data0 = np.loadtxt(filein0)
            xs0 = data0[:, 0]
            ys0 = data0[:, 1]
            x_tr0 = xs0[0]
            y_tr0 = ys0[0]
            data1 = np.loadtxt(filein1)
            xs1 = data1[:, 0]
            ys1 = data1[:, 1]
            x_tr1 = xs1[0]
            y_tr1 = ys1[0]
            assert(y_tr0 < 1e3 and y_tr1 < 1e3)  # at the center of the slab
            V_tr0 = (x_tr1 - x_tr0) / dt
            V_trs.append(V_tr0 / cm)
            ts.append(_time / Myr)
        ax.plot(ts, V_trs, label="trench velocity (cm/yr)")
        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Velocity (cm/yr)')
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
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=1.0e6,
                        help='Time interval, affecting the time steps to visualize')
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
    elif _commend == 'morph_case':
        # slab morphology for a case
        SlabMorphologyCase(arg.inputs, rewrite=1, time_interval=arg.time_interval)
    elif _commend == 'plot_trench':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotTrenchPosition(arg.inputs, time_interval=arg.time_interval)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()