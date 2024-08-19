# -*- coding: utf-8 -*-
r"""(one line description)

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage: plot case running results

        Lib_FOO0_PlotCase plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20

descriptions
    First replate FOO with the name of the project
"""
import numpy as np
import sys, os, argparse, re
import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.ParsePrm import ReadPrmFile
from shilofue.PlotVisit import RunScripts
import shilofue.PlotCombine as PlotCombine
from shilofue.ThDSubduction0.PlotVisit import VISIT_OPTIONS, PREPARE_RESULT_OPTIONS
import shilofue.PlotCase as PlotCase
import shilofue.PlotRunTime as PlotRunTime
import shilofue.PlotStatistics as PlotStatistics
import shilofue.ThermalModel as ThermalModel
import shilofue.Rheology as Rheology
 

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
This scripts generate plots for a single case in TwoDSubduction project\n\
\n\
Examples of usage: \n\
\n\
  - default usage: plot case running results, -t option deals with a time range, default is a whole range.\
-ti option deals with a time interval and is useful for animation.\n\
\n\
        Lib_ThDSubduction0_PlotCase plot_case -i  `pwd`\
 -t 0.0 -t1 0.5e6 -ti 0.1e6\n\
\n\
  - plot cases in a directory (loop), same options as before:\n\
        Lib_ThDSubduction0_PlotCase  plot_case_in_dir -i `pwd`\n\
\n\
  - prepare result of a single step by combining figures:\n\
        Lib_ThDSubduction0_PlotCase prepare_result_step -i `pwd`\n\
\n\
  - plot and then prepare result of a single step by combining figures:\n\
        Lib_ThDSubduction0_PlotCase plot_prepare_result_step -i `pwd`\n\
\n\
  - generate animation: \n\
        Lib_ThDSubduction0_PlotCase animate_case -i `pwd`\n\
\n\
  - generate animation for cases in a directory: \n\
        Lib_ThDSubduction0_PlotCase animate_case_in_dir -i `pwd`\n\
\n\
        ")


def PlotCaseRun(case_path, **kwargs):
    '''
    Plot case run result
    Inputs:
        case_path(str): path to the case
        kwargs:
            time_range
            step(int): if this is given as an int, only plot this step
    Returns:
        -
    '''
    print("PlotCaseRun in ThDSubduction: operating")
    time_interval = kwargs.get('time_interval', None)
    step = kwargs.get('step', None)
    last_step = kwargs.get('last_step', 3)

    # steps to plot: here I use the keys in kwargs to allow different
    # options: by steps, a single step, or the last step
    if type(step) == int:
        kwargs["steps"] = [step]
    elif type(step) == list:
        kwargs["steps"] = step
    elif type(step) == str:
        kwargs["steps"] = step
    else:
        kwargs["last_step"] = last_step

    # get case parameters
    prm_path = os.path.join(case_path, 'output', 'original.prm')
    # plot with paraview
    # initiate class object
    Paraview_Options = VISIT_OPTIONS(case_path)
    # call function
    Paraview_Options.Interpret(**kwargs)
    # ofile = os.path.join('visit_scripts', 'slab_sph.py')
    ofile_list = ['slab.py']
    for ofile_base in ofile_list:
        ofile = os.path.join(case_path, 'paraview_scripts', ofile_base)
        paraview_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts',"ThDSubduction", ofile_base)
        paraview_base_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', 'base.py')  # base.py : base file
        Paraview_Options.read_contents(paraview_base_script, paraview_script)  # this part combines two scripts
        Paraview_Options.substitute()  # substitute keys in these combined file with values determined by Interpret() function
        ofile_path = Paraview_Options.save(ofile, relative=False)  # save the altered script
        print("\t File generated: %s" % ofile_path)
    
    # return the Visit_Options for later usage
    return Paraview_Options


class PLOTTER(PlotCase.PLOTTER):
    '''
    A class for preparing results
    '''
    def __init__(self, module, PlotCaseRuns, **kwargs):
        '''
        Initiation
        module (class) - class to use for generating results at each step
        kwargs(dict):

        '''
        PlotCase.PLOTTER.__init__(self, module, PlotCaseRuns)   


def PrScriptToUse(case_path, default):
    '''
    Return the script to use for preparing results
    This is mainly used for generating figures for animation
    Inputs:
    '''
    pr_script = default
    return pr_script

# todo_r10
def ExportStrenghProfile(case_dir, sp_age, slab_str, strain_rate, o_path=None):
    '''
    export the strength profile from my 3d case
    Inputs:
        case_dir - case directory
        sp_age - age of the subducting plate
        slab_str - maximum yield strength used in model
        strain_rate - strain rate for the analysis
    '''
    # load aspect rheology
    rheology_aspect_json = "/mnt/lochy/ASPECT_DATA/ThDSubduction/EBA_2d_consistent_8_6/eba3d_width61_c22_AR4/configurations/mantle_profile_aspect_v1_HK03_wet_mod_dEdiff0.0000e+00_dEdisl0.0000e+00_dVdiff0.000000e+00_dVdisl0.0000e+00_dAdiff1.0000e+00_dAdisl1.0000e+00.json"
    assert(os.path.isfile(rheology_aspect_json))
    with open(rheology_aspect_json, 'r') as fin:
        rheology_aspect = json.load(fin)
    
    # define global parameters
    g = 9.81
    rho_m = 3300.0 # mantle density
    year = 365.25 * 24 * 3600.0  # s in year
    Myr = 1e6 * year
    
    # pressure profile
    zs_test = np.linspace(0.0, 100e3, 1000)
    Ps = zs_test * rho_m * g
    
    # thermal model
    age = sp_age * Myr
    PlateModel = ThermalModel.PLATE_MODEL(150e3, 1e-6, 273.15, 1673.0)
    Ts = PlateModel.T(zs_test, age)
    
    # diffusion & dislocation creep
    # here I took the parameterization for my 2d cases
    da_file = os.path.join(case_dir, 'output', 'depth_average.txt')
    assert(os.path.isfile(da_file))
    Operator = Rheology.RHEOLOGY_OPR()
    Operator.ReadProfile(da_file)
    
    # compute the composite rheology
    diffusion_aspect = rheology_aspect['diffusion_creep']
    eta_diffusion = Rheology.CreepRheologyInAspectViscoPlastic(diffusion_aspect, strain_rate, Ps, Ts)
    taus_diffusion = 2 * strain_rate * eta_diffusion
    dislocation_aspect = rheology_aspect['dislocation_creep']
    eta_dislocation = Rheology.CreepRheologyInAspectViscoPlastic(dislocation_aspect, strain_rate, Ps, Ts)
    taus_dislocation = 2 * strain_rate * eta_dislocation
    eta_dfds = 1.0 / (1.0 / eta_diffusion + 1.0 / eta_dislocation)
    taus_dfds = 2 * strain_rate * eta_dfds
    
    # brittle yielding
    tau_0 = 50e6
    tau_m = slab_str
    friction = 0.47
    taus_test_brittle = Rheology.CoulumbYielding(Ps, tau_0, friction)
    mask_limit = (taus_test_brittle > tau_m)
    taus_test_brittle[mask_limit] = tau_m
    eta_brittle = taus_test_brittle / 2.0 / strain_rate
    
    # peierls creep
    dV = dislocation_aspect['V']
    creep = Rheology.GetPeierlsRheology("MK10")
    taus_peierls = np.zeros(zs_test.size) # Mpa
    for i in range(zs_test.size):
        taus_peierls[i] = Rheology.PeierlsCreepStress(creep, strain_rate, Ps[i], Ts[i], dV=dV)
    eta_peierls = taus_peierls * 1e6 / 2.0 / strain_rate
    
    # get the stress from the composite rheology
    eta = 1.0 / (1.0/eta_brittle + 1.0/eta_peierls + 1.0/eta_dfds)
    eta1 = np.minimum(eta_peierls, eta_dfds)  # minimum
    eta1 = np.minimum(eta_brittle, eta1)
    eta2 = np.minimum(eta_brittle, 1.0 / (1.0 / eta_peierls + 1.0/eta_dfds)) # DP yield
    eta3 = np.minimum(eta_brittle, 1.0/(1.0 / np.minimum(eta_peierls, eta_dislocation) + 1.0 / eta_diffusion)) # competing Peierls and Dislocation
    eta_nopc = np.minimum(eta_brittle, eta_dfds)
    taus = 2 * strain_rate * eta
    taus1 = 2 * strain_rate * eta1
    taus2 = 2 * strain_rate * eta2
    taus3 = 2 * strain_rate * eta3
    taus_nopc = 2 * strain_rate * eta_nopc
    
    # get the minimum
    taus_minimum = np.minimum(taus_test_brittle, taus_peierls*1e6)
    taus_minimum = np.minimum(taus_minimum, taus_dfds)
    
    # get the minimum
    taus_minimum = np.minimum(taus_test_brittle, taus_peierls*1e6)
    taus_minimum = np.minimum(taus_minimum, taus_dfds)

    if o_path is not None: 
        # plot the depth profile
        fig = plt.figure(tight_layout=True, figsize=(12, 10))
        gs = gridspec.GridSpec(2, 2)
        
        # shear stress vs depth
        ax = fig.add_subplot(gs[0, 0])
        ax.plot(taus2/1e6, zs_test/1e3, label="Composite, DP yield")
        ax.plot(taus_nopc/1e6, zs_test/1e3, label="Composite (no peierls)")
        ax.plot(taus_test_brittle/1e6, zs_test/1e3, 'b--', label="Brittle")
        ax.plot(taus_peierls, zs_test/1e3, 'c--', label="Peierls") # peierls is MPa
        ax.plot(taus_dfds/1e6, zs_test/1e3, 'g--', label="Diff-Disl")
        ax.set_xlabel("Shear Stress (MPa)")
        ax.set_xlim([0, tau_m/1e6 * 1.1])
        ax.invert_yaxis()
        ax.set_ylabel("Depth (km)")
        ax.grid()
        ax.legend()
        
        # viscosity vs depth
        ax = fig.add_subplot(gs[0, 1])
        ax.semilogx(eta2, zs_test/1e3, label="Composite, DP yield")
        ax.semilogx(eta_brittle, zs_test/1e3, 'b--', label="Brittle")
        ax.semilogx(eta_peierls, zs_test/1e3, 'c--', label="Peierls")
        ax.semilogx(eta_dfds, zs_test/1e3, 'g--', label="Diff-Disl")
        ax.set_xlabel("Viscosity (Pa s)")
        ax.set_xlim([np.amin(eta)/10.0, np.amax(eta)*10.0])
        ax.invert_yaxis()
        ax.set_ylabel("Depth (km)")
        ax.grid()
        ax.set_title("strain rate = %.2e" % strain_rate)
        
        # shear stress vs depth
        ax = fig.add_subplot(gs[1, 0])
        ax.plot(taus2/1e6, zs_test/1e3, label="Composite, DP yield")
        ax.plot(taus/1e6, zs_test/1e3, label="Composite, SL yield")
        ax.plot(taus_nopc/1e6, zs_test/1e3, label="Composite (no peierls)")
        ax.plot(taus1/1e6, zs_test/1e3, label="Minimization")
        ax.plot(taus3/1e6, zs_test/1e3, label="Min(Peierls, dislocation)")
        ax.set_xlabel("Shear Stress (MPa)")
        ax.set_xlim([0, tau_m/1e6 * 1.1])
        ax.invert_yaxis()
        ax.set_ylabel("Depth (km)")
        ax.grid()
        ax.legend()

        fig.savefig(o_path)
        print("Saved figure %s" % o_path)
    
    return zs_test, taus2, eta2


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
    parser.add_argument('-t', '--time', type=float,
                        default=None,
                        help='Time')
    parser.add_argument('-t1', '--time1', type=float,
                        default=None,
                        help='Time1')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='step')
    parser.add_argument('-r', '--rewrite', type=int,
                        default=0,
                        help='If rewrite previous result')
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=None,
                        help='Time interval, affecting the time steps to visualize')
    
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    default = os.path.join(ASPECT_LAB_DIR, "files", "ThDSubduction", "figure_step_template.json")

    # commands
    # rearrange the time_range entry
    if arg.time == None or arg.time1 == None:
        time_range = None
    else:
        assert(type(arg.time) == float and type(arg.time1) == float)
        time_range = [arg.time, arg.time1]
    if _commend == 'plot_case':
        PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range, time_interval=arg.time_interval)
        # PlotCase.PlotCaseCombined([PlotCaseRun], arg.inputs, time_range=time_range, time_interval=arg.time_interval)
    elif _commend == 'plot_case_in_dir':
        PlotCase.PlotCaseCombinedDir([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range)
        pass
    elif _commend == 'prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'plot_prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PlotPrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'animate_case':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCaseResults(Plotter.PlotPrepareResultStep, arg.inputs, pr_script, time_interval=arg.time_interval)
    elif _commend == 'animate_case_in_dir':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCombinedDir(Plotter.PlotPrepareResultStep, arg.inputs, pr_script)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()
    else:
        raise ValueError("Invalid command: use -h for help information")

# run script
if __name__ == '__main__':
    main()
