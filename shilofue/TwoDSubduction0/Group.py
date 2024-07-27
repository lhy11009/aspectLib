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
import shilofue.Group as GroupP
from shilofue.Group import CreateGroup, ReadBasicInfoGroup
from shilofue.TwoDSubduction0.Cases import CASE, CASE_OPT
from shilofue.TwoDSubduction0.VtkPp import SLABPLOT
from shilofue.Plot import LINEARPLOT
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec

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
  - default usage, create group: \n\
        Lib_TwoDSubduction0_Group create_group -j ~/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/TwoDSubduction/test_group/test.json\n\
  - Plot diagram \n\
        the -i option is the directory to plot. An additional -o could be append for the output\n\
        note an existing output file will be first parsed in \n\
        python -m shilofue.TwoDSubduction0.Group plot_group_diagram -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT19 -o /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT18/case_summary.txt\n\
  - generate the script for slab morphology\n\
        python -m shilofue.TwoDSubduction0.Group generate_morph_script -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT18/ \n\
\n\
"
        )

def ShowJsonOption():
    Group_Opt = GROUP_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Group_Opt.document_str()
        )


class GROUP_OPT(GroupP.GROUP_OPT):
    pass

class GROUP(GroupP.GROUP):
    pass


class CASE_SUMMARY(GroupP.CASE_SUMMARY):
    '''
    Attributes:
        cases: name of cases
        steps: end steps of cases
        times: end times of cases
        wallclocks: running time of cases on the wall clock
        ab_paths: absolution_paths of cases
        t660s: time the slab tip reaches 660 km
    '''
    def __init__(self, **kwargs):
        '''
        initiation
        Inputs:
            kwargs
        '''
        GroupP.CASE_SUMMARY.__init__(self, **kwargs)
        self.t660s = []
        self.attrs.append("t660s")
        self.sz_methods = []
        self.attrs.append('sz_methods')
        self.sz_thicks = []
        self.attrs.append("sz_thicks")
        self.sz_depths = []
        self.attrs.append("sz_depths")
        self.sz_viscs = []
        self.attrs.append("sz_viscs")
        self.slab_strs = []
        self.attrs.append("slab_strs")
        self.sd_modes = []
        self.attrs.append("sd_modes")
        self.V_sink_avgs = []
        self.attrs.append("V_sink_avgs")
        self.V_plate_avgs = []
        self.attrs.append("V_plate_avgs")
        self.V_ov_plate_avgs = []
        self.attrs.append("V_ov_plate_avgs")
        self.V_trench_avgs = []
        self.attrs.append("V_trench_avgs")
        self.sp_ages = []
        self.attrs.append("sp_ages")
        self.ov_ages = []
        self.attrs.append("ov_ages")

        self.attrs_to_output += self.attrs
        # ['t660s', "sz_thicks", "sz_depths", "sz_viscs", "slab_strs", "sd_modes", "V_sink_avgs", "V_plate_avgs", "V_ov_plate_avgs", "V_trench_avgs", "sp_ages", "ov_ages"]
        # self.headers += ['t660s (yr)', "sz_thicks (m)", "sz_depths (m)", "sz_viscs (Pa s)", "slab_strs (Pa)", "sd_modes", "V_sink_avgs (m/s)", "V_plate_avgs (m/s)", "V_ov_plate_avgs (m/s)", "V_trench_avgs (m/s)", "sp_ages (yr)", "op_ages (yr)"]
    
    def Update(self, **kwargs):
        '''
        Update on properties
        Inputs:
            kwargs:
                actions (list): actions to take
        '''
        actions = kwargs.get('actions', [])

        GroupP.CASE_SUMMARY.Update(self, **kwargs)

        if "t660" in actions:
            # initiate these field
            self.t660s = [-1 for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_t660(i)

        if "shear_zone" in actions:
            # initiate these field
            self.sz_thicks = [-1 for i in range(self.n_case)]
            self.sz_depths = [-1 for i in range(self.n_case)]
            self.sz_viscs = [-1 for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_shear_zone(i)

        if "strength" in actions:
            self.slab_strs = [-1 for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_slab_strength(i)

        if "sz_method" in actions:
            self.sz_methods = [-1 for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_sz_methods(i)

        if "sd_modes" in actions:
            # assign an default value of 1 to the subducting modes
            if len(self.sd_modes) < self.n_case:
                self.sd_modes += [1 for i in range(self.n_case - len(self.sd_modes))]

        if "Vage" in actions:
            # calculate the average velocities
            # initiate these field
            self.V_sink_avgs = [-1 for i in range(self.n_case)]
            self.V_plate_avgs = [-1 for i in range(self.n_case)]
            self.V_ov_plate_avgs = [-1 for i in range(self.n_case)]
            self.V_trench_avgs = [-1 for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_Vavg(i, **kwargs)
        
        if "ages" in actions:
            # ov and sp ages
            self.sp_ages = [-1 for i in range(self.n_case)]
            self.ov_ages = [-1 for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_plate_ages(i)

    # todo_diagram
    def generate_py_scripts(self):
        '''
        generate py script to generate missing slab_morph.txt files
        '''
        py_temp_file = os.path.join(ASPECT_LAB_DIR, 'py_temp.sh')
        py_commands = []
        for i in range(self.n_case):
            if float(self.includes[i]) > 0.0 and self.t660s[i] < 0.0:
                time_interval = 1e5
                case_dir = self.ab_paths[i]
                filename = "slab_morph_t%.2e.txt" % time_interval
                file_path = os.path.join(case_dir, 'vtk_outputs', filename)
                if not os.path.isfile(file_path):
                    py_command = "python -m shilofue.TwoDSubduction0.VtkPp morph_case_parallel -i %s -ti %.2e\n" % (case_dir, time_interval)
                    py_commands.append(py_command)

        if len(py_commands) > 0: 
            with open(py_temp_file, 'w') as fout:
                fout.write("#!/bin/bash\n")
                for py_command in py_commands:
                    fout.write(py_command)
            print("%d cases missing slab morphology outputs, generate morph script: %s" % (len(py_commands), py_temp_file))

    def import_directory(self, _dir, **kwargs):
        '''
        Import from a directory, look for groups and cases
        Inputs:
            _dir (str): directory to import
            kwargs
        '''
        assert(os.path.isdir(_dir))
        GroupP.CASE_SUMMARY.import_directory(self, _dir, **kwargs)

    def update_t660(self, i):
        '''
        update t660
        '''
        case_dir = self.ab_paths[i]
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)
        # use the SLABPLOT class to read the slab_morph.txt file
        # and get the t660
        SlabPlot = SLABPLOT('foo')
        try:
            t660 = SlabPlot.GetTimeDepthTip(case_dir, 660e3, filename="slab_morph_t1.00e+05.txt")
        except SLABPLOT.SlabMorphFileNotExistError:
            t660 = -1.0
        self.t660s[i] = t660
    
    def update_Vavg(self, i, **kwargs):
        '''
        update the average velocities
        '''
        case_dir = self.ab_paths[i]
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)
        t1_method = kwargs.get("t1_method", "value")
        assert(t1_method in ["value", 't660'])
        # t1 could either be assigned by value or by a factor * t660
        if t1_method == "value":
            # by default, a 30 Ma time window is used to average the velocity
            t1 = kwargs.get("t1", 30e6)
        elif t1_method == "t660":
            # by default, 5.0 * t660 is the time window used to average the velocity
            t1_factor = kwargs.get("t1_factor", 5.0)
            t660 = self.t660s[i]
            if t660 > 0.0:
                t1 = 5.0 * t660
            else:
                # in case the t660 is not computed yet, take the -1.0 as the default
                # this will lead to an invalid time range and trigger the value -2.0 in
                # calculating the velocities.
                t1 = -1.0
        # use the SLABPLOT class to read the slab_morph.txt file
        # and get the t660
        SlabPlot = SLABPLOT('foo')
        try:
            t660 = SlabPlot.GetTimeDepthTip(case_dir, 660e3, filename="slab_morph_t1.00e+05.txt")
            # sink velocity, sp plate velocity and the trench velocity averaged in a time range
            V_sink_avg, V_plate_avg, V_ov_plate_avg, V_trench_avg = SlabPlot.GetAverageVelocities(case_dir, t660, t1, filename="slab_morph_t1.00e+05.txt")
        except SLABPLOT.SlabMorphFileNotExistError:
            print("morph file not exsit in dir: %s" % case_dir) # debug
            V_sink_avg = -1.0
            V_plate_avg = -1.0
            V_ov_plate_avg = -1.0
            V_trench_avg = -1.0
        self.V_sink_avgs[i] = V_sink_avg
        self.V_plate_avgs[i] = V_plate_avg
        self.V_ov_plate_avgs[i] = V_ov_plate_avg
        self.V_trench_avgs[i] = V_trench_avg
        
    def update_shear_zone(self, i):
        '''
        Update shear zone properties
        ''' 
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.sz_viscs[i] = Visit_Options.options["SHEAR_ZONE_CONSTANT_VISCOSITY"]
            self.sz_depths[i] = Visit_Options.options["SHEAR_ZONE_CUTOFF_DEPTH"]
            self.sz_thicks[i] = Visit_Options.options["INITIAL_SHEAR_ZONE_THICKNESS"]
        except FileNotFoundError:
            self.sz_viscs[i] = -1.0
            self.sz_depths[i] = -1.0
            self.sz_thicks[i] = -1.0

    def update_slab_strength(self, i):
        '''
        Update slab strength properties
        ''' 
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.slab_strs[i] = Visit_Options.options["MAXIMUM_YIELD_STRESS"]
        except FileNotFoundError:
            self.slab_strs[i] = -1.0
    
    def update_sz_methods(self, i):
        '''
        Update method of shear zone viscosity
        '''
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.sz_methods[i] = Visit_Options.options["SHEAR_ZONE_METHOD"]
        except FileNotFoundError:
            self.sz_methods[i] = -1
        pass

    def update_plate_ages(self, i):
        '''
        update plate ages
        '''
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.sp_ages[i] = Visit_Options.options["SP_AGE"]
            self.ov_ages[i] = Visit_Options.options["OV_AGE"]
        except FileNotFoundError:
            self.sp_ages[i] = -1.0
            self.ov_ages[i] = -1.0

    def plot_diagram_t660(self, **kwargs):
        '''
        plot a diagram of t660
        '''
        fig_path = kwargs.get('fig_path', None)

        fig, ax = plt.subplots()
        x_name = "sz_thicks"
        y_name = "sz_depths"
        z_name = "t660s"

        # mesh data
        includes = np.array([int(i) for i in getattr(self, "includes")])
        SZTT = np.array(getattr(self, x_name))
        SZDD= np.array(getattr(self, y_name))
        tt660 = np.array(getattr(self, z_name))
        sd_modes = np.array([int(i) for i in getattr(self, "sd_modes")])

        # include only case with subducting mode 1
        mask = (includes > 0) & (sd_modes == 1)

        # generate the plot
        # h = ax.pcolormesh(SZTT/1e3, SZDD/1e3, tt660)
        h = ax.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, c=tt660[mask]/1e6, vmin=0.0, vmax=5.0)
        ax.set_xlabel("shear zone thickness (km)")
        ax.set_ylabel("shear zone depth (km)")
        fig.colorbar(h, ax=ax, label='t_660 (Ma)')

        # save figure
        if fig_path is not None:
            fig.savefig(fig_path)
            print("figure generated: ", fig_path)

    def plot_velocities_tr(self, **kwargs):
        '''
        plot a diagram of velocities
        '''
        fig_path = kwargs.get('fig_path', None)

        x_name = "sz_thicks"
        y_name = "sz_depths"

        # mesh data
        includes = np.array([int(i) for i in getattr(self, "includes")])
        SZTT = np.array(getattr(self, x_name))
        SZDD= np.array(getattr(self, y_name))
        VVTR = np.array(getattr(self, "V_trench_avgs"))
        VVSP = np.array(getattr(self, "V_plate_avgs"))
        VVSK = np.array(getattr(self, "V_sink_avgs"))
        VV = VVSP - VVTR
        sd_modes = np.array([int(i) for i in getattr(self, "sd_modes")])

        # include only case with subducting mode 1
        mask1 = (includes > 0) & (sd_modes == 1)
        mask2 = (VVTR > -1.0)
        mask = mask1&mask2

        # generate axis
        fig = plt.figure(tight_layout=True, figsize=(12, 10))
        gs = gridspec.GridSpec(2,2)
        area_scale = 200.0 # scaling for the area of scatter points

        # generate the plot
        # figure 1: plate velocity 
        ax0 = fig.add_subplot(gs[0, 0])
        area = (area_scale*VVSP[mask])**2.0
        h = ax0.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area, c=100.0*VVSP[mask])
        ax0.set_xlabel("shear zone thickness (km)")
        ax0.set_ylabel("shear zone depth (km)")
        ax0.set_title("V_plate_avgs")
        fig.colorbar(h, ax=ax0, label='V (cm/yr)')

        # figure 2: trench velocity
        # negative value converted to positive here
        # area = (200.0*VV[mask])**2.0
        # h = ax.scatter(SZTT[mask&mask1]/1e3, SZDD[mask&mask1]/1e3, s=area, c=100.0*VV[mask&mask1], vmin=0.0, vmax=5.0)
        ax1 = fig.add_subplot(gs[0, 1])
        area = (area_scale*VVTR[mask])**2.0
        h = ax1.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area, c=-100.0*VVTR[mask], vmin=0.0, vmax=5.0)
        # h = ax.scatter(SZTT[mask&mask1]/1e3, SZDD[mask&mask1]/1e3, c=-100.0*VVTR[mask&mask1], vmin=0.0, vmax=5.0)
        ax1.set_xlabel("shear zone thickness (km)")
        ax1.set_ylabel("shear zone depth (km)")
        ax1.set_title("V_trench_avgs")
        fig.colorbar(h, ax=ax1, label='V (cm/yr)')

        # figure 3: trench and plate velocity
        ax2 = fig.add_subplot(gs[1, 0])
        area = (area_scale*VVTR[mask])**2.0
        area1 = (area_scale*VV[mask])**2.0
        ax2.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area1, c='w', edgecolors='k') # plot the bigger points of total velocity
        h = ax2.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area, c=-100.0*VVTR[mask], vmin=0.0, vmax=5.0) # plot the smaller points of trench velocity
        ax2.set_xlabel("shear zone thickness (km)")
        ax2.set_ylabel("shear zone depth (km)")
        ax2.set_title("V_trench_avgs")
        fig.colorbar(h, ax=ax2, label='V (cm/yr)')
        
        # figure 4: trench and plate velocity - annotation
        n_a = 5
        dv = 0.02
        ax3 = fig.add_subplot(gs[1, 1])
        velocities = np.linspace(dv, dv*n_a, n_a)
        xs = np.ones(n_a)
        ys = np.linspace(1.0, 1.0*n_a, n_a)
        area = (area_scale*velocities)**2.0
        ax3.scatter(xs, ys, s=area, c='w', edgecolors='k')
        for i in range(n_a):
            ax3.annotate( "%.1f cm/yr" % (100.0*dv*(i+1)), xy=(1.0, 1.0*(i+1)),
                    xytext=(0, 12),  # 4 points vertical offset.
                    textcoords='offset points',
                    ha='center', va='bottom')
        ax3.set_xlim([0.0, 6.0])
        ax3.set_ylim([0.5, n_a + 0.5])
        ax3.set_axis_off()

        # save figure
        if fig_path is not None:
            fig.savefig(fig_path)
            print("figure generated: ", fig_path)

    def plot_velocities_tr_publication(self, fig_path):
        '''
        plot a diagram of velocities
        '''
        x_name = "sz_thicks"
        y_name = "sz_depths"

        # mesh data
        includes = np.array([int(i) for i in getattr(self, "includes")])
        SZTT = np.array(getattr(self, x_name))
        SZDD= np.array(getattr(self, y_name))
        VVTR = np.array(getattr(self, "V_trench_avgs"))
        VVSP = np.array(getattr(self, "V_plate_avgs"))
        VVSK = np.array(getattr(self, "V_sink_avgs"))
        VV = VVSP - VVTR
        sd_modes = np.array([int(i) for i in getattr(self, "sd_modes")])

        # include only case with subducting mode 1
        mask1 = (includes > 0) & (sd_modes == 1)
        mask2 = (VVTR > -1.0)
        mask = mask1&mask2

        # generate axis
        fig = plt.figure(tight_layout=True, figsize=(18, 7.5))
        gs = gridspec.GridSpec(1, 2)
        area_scale = 200.0 # scaling for the area of scatter points

        # figure: trench and plate velocity
        ax0 = fig.add_subplot(gs[0, 0])
        area = (area_scale*VVTR[mask])**2.0
        area1 = (area_scale*VV[mask])**2.0
        ax0.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area1, c='w', edgecolors='k') # plot the bigger points of total velocity
        h = ax0.scatter(SZTT[mask]/1e3, SZDD[mask]/1e3, s=area, c=-100.0*VVTR[mask], vmin=0.0, vmax=5.0) # plot the smaller points of trench velocity
        ax0.set_xlabel("shear zone thickness (km)")
        ax0.set_ylabel("shear zone depth (km)")
        fig.colorbar(h, ax=ax0, label='V (cm/yr)')
        
        # figure annotation
        n_a = 5
        dv = 0.02
        ax1 = fig.add_subplot(gs[0, 1])
        velocities = np.linspace(dv, dv*n_a, n_a)
        xs = np.ones(n_a)
        ys = np.linspace(1.0, 1.0*n_a, n_a)
        area = (area_scale*velocities)**2.0
        ax1.scatter(xs, ys, s=area, c='w', edgecolors='k')
        for i in range(n_a):
            ax1.annotate( "%.1f cm/yr" % (100.0*dv*(i+1)), xy=(1.0, 1.0*(i+1)),
                    xytext=(0, 12),  # 4 points vertical offset.
                    textcoords='offset points',
                    ha='center', va='bottom')
        ax1.set_xlim([0.0, 6.0])
        ax1.set_ylim([0.5, n_a + 0.5])
        ax1.set_axis_off()
            
        fig.savefig(fig_path)
        print("figure generated: ", fig_path)

    def plot_velocities_ages_publication(self, fig_path):
        '''
        plot a diagram of velocities
        '''
        x_name = "sp_ages"
        y_name = "ov_ages"

        # mesh data
        includes = np.array([int(i) for i in getattr(self, "includes")])
        SPAG = self.export(x_name)
        OVAG = self.export(y_name)
        VVTR = self.export("V_trench_avgs")
        VVSP = self.export("V_plate_avgs")
        VVSK = self.export("V_sink_avgs")
        VV = VVSP - VVTR
        sd_modes = np.array([int(i) for i in getattr(self, "sd_modes")])

        # include only case with subducting mode 1
        mask1 = (includes > 0) & (sd_modes == 1)
        mask2 = (VVTR > -1.0)
        mask = mask1&mask2

        # generate axis
        fig = plt.figure(tight_layout=True, figsize=(18, 7.5))
        gs = gridspec.GridSpec(1, 2)
        area_scale = 200.0 # scaling for the area of scatter points

        # figure: trench and plate velocity
        ax0 = fig.add_subplot(gs[0, 0])
        area = (area_scale*VVTR[mask])**2.0
        area1 = (area_scale*VV[mask])**2.0
        ax0.scatter(SPAG[mask]/1e6, OVAG[mask]/1e6, s=area1, c='w', edgecolors='k') # plot the bigger points of total velocity
        h = ax0.scatter(SPAG[mask]/1e6, OVAG[mask]/1e6, s=area, c=-100.0*VVTR[mask], vmin=0.0, vmax=5.0) # plot the smaller points of trench velocity
        ax0.set_xlabel("sp age (Ma)")
        ax0.set_ylabel("ov age (Ma)")
        fig.colorbar(h, ax=ax0, label='V (cm/yr)')
        
        # figure annotation
        n_a = 5
        dv = 0.02
        ax1 = fig.add_subplot(gs[0, 1])
        velocities = np.linspace(dv, dv*n_a, n_a)
        xs = np.ones(n_a)
        ys = np.linspace(1.0, 1.0*n_a, n_a)
        area = (area_scale*velocities)**2.0
        ax1.scatter(xs, ys, s=area, c='w', edgecolors='k')
        for i in range(n_a):
            ax1.annotate( "%.1f cm/yr" % (100.0*dv*(i+1)), xy=(1.0, 1.0*(i+1)),
                    xytext=(0, 12),  # 4 points vertical offset.
                    textcoords='offset points',
                    ha='center', va='bottom')
        ax1.set_xlim([0.0, 6.0])
        ax1.set_ylim([0.5, n_a + 0.5])
        ax1.set_axis_off()
            
        fig.savefig(fig_path)
        print("figure generated: ", fig_path)
    
    def plot_velocities_sink(self, **kwargs):
        '''
        plot a diagram of sink velocities
        '''
        fig_path = kwargs.get('fig_path', None)

        fig, ax = plt.subplots()
        x_name = "sz_thicks"
        y_name = "sz_depths"
        z_name = kwargs.get("z_name", "V_sink_avgs")

        # mesh data
        includes = np.array([int(i) for i in getattr(self, "includes")])
        SZTT = np.array(getattr(self, x_name))
        SZDD= np.array(getattr(self, y_name))
        VVSK = np.array(getattr(self, z_name))
        sd_modes = np.array([int(i) for i in getattr(self, "sd_modes")])

        # include only case with subducting mode 1
        mask = (includes > 0) & (sd_modes == 1)
        mask1 = (VVSK > -1.0)

        # generate the plot
        # negative value converted to positive here
        h = ax.scatter(SZTT[mask&mask1]/1e3, SZDD[mask&mask1]/1e3, c=100.0*VVSK[mask&mask1])
        ax.set_xlabel("shear zone thickness (km)")
        ax.set_ylabel("shear zone depth (km)")
        ax.set_title(z_name)
        fig.colorbar(h, ax=ax, label='V (cm/yr)')

        # save figure
        if fig_path is not None:
            fig.savefig(fig_path)
            print("figure generated: ", fig_path)


def PlotGroupDiagram(group_dir, **kwargs):
    '''
    pass
    '''
    assert(os.path.isdir(group_dir))
    
    o_path = kwargs.get('o_path', "")
    if o_path == "":
        o_path = os.path.join(group_dir, 'case_summary.txt')
    
    fig_t660_path = os.path.join(group_dir, 'img', 't660s.png')
    if not os.path.isdir(os.path.dirname(fig_t660_path)):
      os.mkdir(os.path.dirname(fig_t660_path))
    
    fig_vtr_path = os.path.join(group_dir, 'img', 'vtr.png')
    
    Case_Summary = CASE_SUMMARY(VISIT_OPTIONS=VISIT_OPTIONS)

    # import old result if it exists 
    if os.path.isfile(o_path):
        Case_Summary.import_txt(o_path)
    
    # import the new directory
    Case_Summary.import_directory(group_dir, actions=['t660', 'shear_zone', "sd_modes", "Vage", "ages"])

    # output file output 
    Case_Summary.write_file(o_path)
    
    # plot the t660 diagram
    Case_Summary.plot_diagram_t660(fig_path=fig_t660_path)

    # plot the velocity diagram: trench
    Case_Summary.plot_velocities_tr(fig_path=fig_vtr_path)


def GenerateMorphScript(group_dir):
    '''
    Generate the script to plot morphology
    Inputs:
        group_dir - directory of cases
    '''
    time_interval = 0.1e6
    
    # read the list of cases
    case_list, _, _, _ = ReadBasicInfoGroup(group_dir)
    
    # write a .sh file to run command in the system
    py_temp_file = os.path.join(ASPECT_LAB_DIR, 'py_temp.sh')
    py_commands = []
    
    for _case in case_list:
      case_dir = os.path.join(group_dir, _case)
      filename = "slab_morph_t%.2e.txt" % time_interval
      file_path = os.path.join(case_dir, 'vtk_outputs', filename)
      if not os.path.isfile(file_path):
        print("file_path: ", file_path) # debug
        # print("filename: ", filename) # debug
        py_command = "python -m shilofue.TwoDSubduction0.VtkPp morph_case_parallel -i %s -ti %.2e\n" % (case_dir, time_interval)
        py_commands.append(py_command)
    
    with open(py_temp_file, 'w') as fout:
      fout.write("#!/bin/bash\n")
      for py_command in py_commands:
        fout.write(py_command)
    
    print("GenerateMorphScript: %s" % py_temp_file)


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
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='A json file')
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
    elif (_commend in ['--json_option', '-jo']):
        # json options
        ShowJsonOption()
    elif _commend == 'create_group':
        CreateGroup(arg.json, CASE, CASE_OPT)
    elif _commend == 'plot_group_diagram':
        PlotGroupDiagram(arg.inputs, o_path=arg.outputs)
    elif _commend == 'generate_morph_script':
        GenerateMorphScript(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()