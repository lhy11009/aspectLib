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
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, VISIT_OPTIONS
import shilofue.Utilities as Utilities

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def Usage():
    print("\
This scripts analyze slab morphology\n\
\n\
Examples of usage: \n\
\n\
  - generate slab_morph.txt: \n\
    (what this does is looping throw all visualizing steps, so it takes time)\n\
    python -m shilofue.TwoDSubduction0.PlotSlab morph_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8 \n\
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
    def PlotSlabSurface(self, fileout):
        fig, ax = plt.subplots()
        ax.plot(self.rs * np.cos(self.thetas), self.rs * np.sin(self.thetas), label='slab surface')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        plt.savefig(fileout)
        print("plot slab surface: ", fileout)
    
    def PlotMorph(self, case_dir):
        # path
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        morph_dir = os.path.join(img_dir, 'morphology')
        if not os.path.isdir(morph_dir):
            os.mkdir(morph_dir)
        # read parameters
        Ro = float(self.prm['Geometry model']['Chunk']['Chunk outer radius'])
        # read data
        slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
        assert(os.path.isfile(slab_morph_file))
        self.ReadHeader(slab_morph_file)
        self.ReadData(slab_morph_file)
        col_pvtu_step = self.header['pvtu_step']['col']
        col_pvtu_time = self.header['time']['col']
        col_pvtu_trench = self.header['trench']['col']
        col_pvtu_slab_depth = self.header['slab_depth']['col']
        pvtu_steps = self.data[:, col_pvtu_step]
        times = self.data[:, col_pvtu_time]
        trenches = self.data[:, col_pvtu_trench]
        trenches_migration_length = (trenches - trenches[0]) * Ro  # length of migration
        slab_depthes = self.data[:, col_pvtu_slab_depth]
        # todo
        trench_velocities = np.gradient(trenches_migration_length, times)
        sink_velocities = np.gradient(slab_depthes, times)
        # trench velocity
        # start figure
        fig = plt.figure(tight_layout=True) 
        fig.subplots_adjust(hspace=0)
        gs = gridspec.GridSpec(2, 2) 
        # 1: trench & slab movement
        ax = fig.add_subplot(gs[0, :]) 
        ax_tx = ax.twinx()
        _color = 'tab:blue'
        lns0 = ax.plot(times/1e6, trenches_migration_length/1e3, '-', color=_color, label='trench position (km)')
        ax.set_ylabel('Trench Position (km)', color=_color)
        ax.tick_params(axis='y', labelcolor=_color)
        ax.grid()
        lns1 = ax_tx.plot(times/1e6, slab_depthes/1e3, 'k-', label='slab depth (km)')
        ax_tx.set_ylabel('Slab Depth (km)')
        lns = lns0 + lns1
        labs = [I.get_label() for I in lns]
        # ax.legend(lns, labs)
        # 2:
        # todo
        ax = fig.add_subplot(gs[1, :]) 
        _color = 'tab:blue'
        ax.plot(times/1e6, 0.0 * np.zeros(times.shape), 'k--')
        lns0 = ax.plot(times/1e6, trench_velocities*1e2, '-', color=_color, label='trench velocity (cm/yr)')
        ax.plot(times/1e6, sink_velocities*1e2, 'k-', label='trench velocity (cm/yr)')
        ax.set_ylim((-10, 20))
        ax.set_xlabel('Times (Myr)')
        ax.set_ylabel('Velocity (cm/yr)')
        ax.grid()
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
    output_dir = os.path.join(case_dir, 'vtk_outputs')
    output_file = os.path.join(output_dir, 'slab_morph.txt')
    vtk_option_path, _time, step = PrepareVTKOptions(case_dir, 'TwoDSubduction_SlabAnalysis', vtk_step=pvtu_step)
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    slab_outputs = slab_morph(_stdout)
    # header
    file_header = "# %s\n# %s\n# %s\n# %s\n# %s\n" % \
    ("1: pvtu_step", "2: step", "3: time (yr)", "4: trench (rad)", "5: slab depth (m)")
    # output string
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e\n" % (pvtu_step, step, _time, slab_outputs['trench_theta'], slab_outputs['slab_depth'])
    # remove old file
    is_new = kwargs.get('new', False)
    if is_new and os.path.isfile(output_file):
        os.remove(output_file)
    # output data
    if not os.path.isfile(output_file):
        with open(output_file, 'w') as fout:
            fout.write(file_header)
            fout.write(outputs)
        print('Create output: %s' % output_file)
    else:
        with open(output_file, 'a') as fout:
            fout.write(outputs)
        print('Update output: %s' % output_file)


def vtk_and_slab_morph_case(case_dir, **kwargs):
    '''
    run vtk and get outputs for every snapshots
    Inputs:
        kwargs:
            rewrite: if rewrite previous results
    '''
    # get all available snapshots
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    available_pvtu_snapshots = Utilities.string2list(Visit_Options.options['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'])
    available_pvtu_steps = [i - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) for i in available_pvtu_snapshots]
    # get where previous session ends
    SlabPlot = SLABPLOT('slab')
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
    if os.access(slab_morph_file, os.R_OK):
        # read previous results if they exist
        SlabPlot.ReadHeader(slab_morph_file)
        SlabPlot.ReadData(slab_morph_file)
        col_pvtu_step = SlabPlot.header['pvtu_step']['col']
        last_pvtu_step = int(SlabPlot.data[-1, col_pvtu_step])
    else:
        last_pvtu_step = -1
    # get slab morphology
    if_rewrite = kwargs.get('rewrite', 0)
    for pvtu_step in available_pvtu_steps:
        if pvtu_step <= last_pvtu_step and if_rewrite == 0:
            # skip existing steps
            continue
        if pvtu_step == 0:
            # start new file with the 0th step
            vtk_and_slab_morph(case_dir, pvtu_step, new=True)
        else:
            vtk_and_slab_morph(case_dir, pvtu_step)
    pass
    

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
    parser.add_argument('-s', '--step', type=str,
                        default='0',
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
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph_case(arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        prm_file = os.path.join(arg.inputs, 'case.prm')
        SlabPlot.ReadPrm(prm_file)
        SlabPlot.PlotMorph(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()