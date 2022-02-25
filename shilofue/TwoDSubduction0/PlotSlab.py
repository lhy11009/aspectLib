# -*- coding: utf-8 -*-
r"""Plot Slab Morphology

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
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec
from shilofue.Plot import LINEARPLOT
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
from joblib import Parallel, delayed
import multiprocessing

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
This scripts analyze slab morphology\n\
\n\
Examples of usage: \n\
\n\
  - plot the contour of slab of a give step and interact with the plot\n\
    python -m shilofue.TwoDSubduction0.PlotSlab plot_morph_contour_step -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT1/eba_cdpt_SA80.0_OA40.0 -s 12\n\
\n\
  - generate slab_morph.txt: \n\
    (what this does is looping throw all visualizing steps, so it takes time)\n\
    the -r option rewrite previous slab_morph.txt, otherwise new results will be appended to that previous file\n\
    python -m shilofue.TwoDSubduction0.PlotSlab morph_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8 \n\
    -r 1 \n\
\n\
    Options: \n\
        -r: 0(default), 1 - rewrite previous slab_morph.txt\n\
\n\
  - plot trench movement: \n\
    (note you have to have a slab_morph.txt generated) \n\
    python -m shilofue.TwoDSubduction0.PlotSlab plot_morph -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8 \n\
\n\
        ")

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

    def PlotSlabSurface(self, fileout):
        fig, ax = plt.subplots()
        ax.plot(self.rs * np.cos(self.thetas), self.rs * np.sin(self.thetas), label='slab surface')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        plt.savefig(fileout)
        print("plot slab surface: ", fileout)
    
    def ReadWedgeT(self, case_dir, min_pvtu_step, max_pvtu_step, **kwargs):
        # todo
        time_interval_for_slab_morphology = 0.5e6  # hard in
        i = 0
        initial_adaptive_refinement = int(self.prm['Mesh refinement']['Initial adaptive refinement'])
        Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
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
                depthes = Ro - rs # compute depth
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
        Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
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
        trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        # trench velocity
        # start figure
        fig = plt.figure(tight_layout=True) 
        fig.subplots_adjust(hspace=0)
        gs = gridspec.GridSpec(3, 3) 
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
        ax.set_ylim((-10, 20))
        ax.set_ylabel('Velocity (cm/yr)')
        ax.grid()
        ax.tick_params(axis='x', labelbottom=False) # labels along the bottom edge are off
        # 3: wedge temperature
        depthes, Ts = self.ReadWedgeT(case_dir, int(pvtu_steps[0]), int(pvtu_steps[-1]))
        tt, dd = np.meshgrid(times, depthes)
        ax = fig.add_subplot(gs[2, 0:2]) 
        h = ax.pcolormesh(tt/1e6,dd/1e3,Ts, shading='gouraud') 
        ax.invert_yaxis()
        ax.set_xlim((times[0]/1e6, times[-1]/1e6))  # set x limit
        ax.set_xlabel('Times (Myr)')
        ax.set_ylabel('Depth (km)')
        ax = fig.add_subplot(gs[2, 2])
        ax.axis('off') 
        fig.colorbar(h, ax=ax, label='T (K)') 
        fig.tight_layout()
        # save figure
        o_path = os.path.join(morph_dir, 'trench.png')
        plt.savefig(o_path)
        print("%s: figure %s generated" % (Utilities.func_name(), o_path))



def slab_morph(inputs):
    '''
    Plot slab morphology
    Inputs:
        inputs (str): the outputs from another executable (slab analysis)
    '''
    outputs = {}
    outputs['trench_theta'] = float(Utilities.re_read_variable_from_string(inputs, 'Trench theta', ':'))
    outputs['slab_depth'] = float(Utilities.re_read_variable_from_string(inputs, 'Slab depth', ':'))
    outputs['theta100'] = float(Utilities.re_read_variable_from_string(inputs, '100km theta', ':'))
    outputs['dip100'] = float(Utilities.re_read_variable_from_string(inputs, '100km dip', ':'))
    return outputs


def vtk_and_slab_morph(case_dir, pvtu_step, **kwargs):
    '''
    run vtk and read in slab morph
    Inputs:
        case_dir (str): case directory
        pvtu_step (int): time step
        kwargs:
            new: remove old output file
    '''
    print("pvtu_step: %s\n" % str(pvtu_step))
    vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    vtk_step=pvtu_step, include_step_in_filename=True, generate_horiz=True)
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    slab_outputs = slab_morph(_stdout)
    # output string
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e%-14.4e\n"\
    % (pvtu_step, step, _time, slab_outputs['trench_theta'], slab_outputs['slab_depth'], slab_outputs['dip100'])
    return pvtu_step, outputs


def vtk_and_slab_morph_case(case_dir, **kwargs):
    '''
    run vtk and get outputs for every snapshots
    Inputs:
        kwargs:
            rewrite: if rewrite previous results
    '''
    # get all available snapshots
    # the interval is choosen so there is no high frequency noises
    time_interval_for_slab_morphology = 0.5e6
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # call get_snaps_for_slab_morphology, this prepare the snaps with a time interval in between.
    available_pvtu_snapshots= Visit_Options.get_snaps_for_slab_morphology(time_interval=time_interval_for_slab_morphology)
    available_pvtu_steps = [i - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) for i in available_pvtu_snapshots]
    # get where previous session ends
    SlabPlot = SLABPLOT('slab')
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
    if_rewrite = kwargs.get('rewrite', 0)
    # Read previous result
    last_pvtu_step = -1
    if os.access(slab_morph_file, os.R_OK) and if_rewrite == 0:
        # read previous results if they exist
        SlabPlot.ReadHeader(slab_morph_file)
        SlabPlot.ReadData(slab_morph_file)
        if SlabPlot.HasData():
            col_pvtu_step = SlabPlot.header['pvtu_step']['col']
            last_pvtu_step = int(SlabPlot.data[-1, col_pvtu_step])
    # Initiation Wrapper class for parallel computation
    ParallelWrapper = PARALLEL_WRAPPER_FOR_VTK('slab_morph', vtk_and_slab_morph, last_pvtu_step=last_pvtu_step, if_rewrite=if_rewrite)
    ParallelWrapper.configure(case_dir)  # assign case directory
    # Remove previous file
    if if_rewrite:
        if os.path.isfile(slab_morph_file):
            print("%s: Delete old slab_morph.txt file." % Utilities.func_name())
            os.remove(slab_morph_file)  # delete slab morph file
        ParallelWrapper.delete_temp_files(available_pvtu_steps)  # delete intermediate file if rewrite
    num_cores = multiprocessing.cpu_count()
    for pvtu_step in available_pvtu_steps:
        ParallelWrapper(pvtu_step)  # debug
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


def plot_morph_contour_step(case_dir, step):
    '''
    plot slab morphology contour
    '''
    file_in_path = os.path.join(case_dir, 'vtk_outputs', 'contour_slab_%05d.txt' % (step))
    assert(os.path.isfile(file_in_path))
    ## load data
    data = np.loadtxt(file_in_path)
    fig, ax = plt.subplots()
    ax.plot(data[:, 0], data[:, 1], 'b.')
    ax.set_xlim([0, 6371e3])
    ax.set_ylim([0, 6371e3])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    temp_dir = os.path.join(case_dir, 'temp_output')
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    file_out = os.path.join(temp_dir, "slab_contour_s%05d.png" % step)
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
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='Timestep')
    parser.add_argument('-r', '--rewrite', type=int,
                        default=0,
                        help='If rewrite previous result')
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
    elif _commend == 'morph_step':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph(arg.inputs, int(arg.step))
    elif _commend == 'plot_morph_contour_step':
        plot_morph_contour_step(arg.inputs, arg.step)
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph_case(arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotMorph(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()