# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - convert parameters of some type of rheology to aspect input:

        python -m shilofue.Rheology convert_to_ASPECT -r HK03 -E 1 -j temp.json

        options for "-r":
            HK03 - Hirth & Kohlstedt 2003
            AB17 - Arredondo & Billen 2017
        
        -j: path to a json file to save the output

  - compute the viscosity by first selecting a rheology

        python -m shilofue.Rheology compute_creep_viscosity -r AB17 -P 10e9 -T 1673 -S 1e-15 -E 1

        -E: use effective strain rate, a variable F will be computed and applied as prefactor

  - compute the viscosity with ASPECT's formulation

        python -m shilofue.Rheology compute_ASPECT_viscosity -j temp.json -P 10e9 -T 1673 -S 1e-15

        -j: parameters of rheology are loaded from a json file

        there is no '-E' option, as F should be incorporated in A in ASPECT's rheology

  - plot along a profile from aspect to check different parameterazition

        python -m shilofue.Rheology plot_along_aspect_profile -r HK03 -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/output/depth_average.txt -im MB

        -im : LHY(default) or MB, use my implementation of formula or Magali's (in script flow_law_function)
  
  - plot along a profile from aspect, using rheology parameterization from a prm file to check implementation in code

        python -m shilofue.Rheology plot_along_aspect_profile_with_json 
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/output/depth_average.txt 
        -j /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/case.prm

        -i: depth_average file
        -j: prm file / json file, in case of a json file, it contains diffusion + dislocation creep(as output from the convert_to_ASPECT command)

descriptions
""" 

import json
import os
import sys
import math
import argparse
import numpy as np
from matplotlib import pyplot as plt
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
from shilofue.flow_law_functions import visc_diff_HK
from shilofue.ParsePrm import ParseFromDealiiInput, UpperMantleRheologyViscoPlastic
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D

R = 8.314

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')

class CheckValueError(Exception):
    pass


class RHEOLOGY_PRM():
    """
    class for rheologies
    """
    def __init__(self):
        '''
        Initiation, initiate rheology parameters
        '''
        # dislocation creep in Hirth & Kohlstedt 2003
        self.HK03_disl = \
            {
                "A": 90,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03_diff = \
            {
                "A" : 1.0e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 335e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # dislocation creep in Arredondo & Billen 2017
        self.AB17_disl = \
            {
                "A": 2.57e-20,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 496e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # diffusion creep in Arredondo & Billen 2017
        self.AB17_diff = \
            {
                "A" : 2.85e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 317e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # modify dislocation creep in Hirth & Kohlstedt 2003
        self.HK03v1_disl = \
            {
                "A": 0.9,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03v1_diff = \
            {
                "A" : 1.0e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 335e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }


class RHEOLOGY_OPR():
    '''
    rheology operation, do some complex staff
    Attributes:
        RheologyPrm: an initiation of the class RHEOLOGY_PRM
        depths(ndarray): depth profile
        pressures: pressure profile
        temperatures: temperature profile
    '''
    pass

    def __init__(self):
        self.RheologyPrm = RHEOLOGY_PRM()
        self.depths = None
        self.pressures = None
        self.tempertures = None
        pass
    
    def ReadProfile(self, file_path):
        self.depths, self.pressures, self.temperatures = ReadAspectProfile(file_path)

    def ConstrainRheology(self, **kwargs):
        '''
        varying around a give rheology with variation with applied constraints
        inputs: 
            kwargs(dict):
                rheology: type of initial rheology
        '''
        constrained_rheologies = []
        constrained_ds = []
        constrained_Vdisls = []
        constrained_Vdiffs = []
        rheology = kwargs.get('rheology', 'HK03')
        diffusion_creep = getattr(self.RheologyPrm, rheology + "_diff")
        dislocation_creep = getattr(self.RheologyPrm, rheology + "_disl")
        strain_rate = 1e-15
        
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # grain size
        d_mean = kwargs.get('d', 0.75e4)

        Vdiff_sigma = 4e-6
        Vdisl_sigma = 11e-6
        d_sigma = 5e3

        # random walk
        eta660range = [4e20, 1e21]
        N = 10000
        Vdiffs = np.random.normal(diffusion_creep['V'], Vdiff_sigma, N)
        Vdisls = np.random.normal(dislocation_creep['V'], Vdisl_sigma, N)
        ds = np.random.normal(d_mean, d_sigma, N)
        for i in range(N):
            Vdiff = Vdiffs[i]
            Vdisl = Vdisls[i]
            d = ds[i]
            if Vdiff < 0.0 or Vdisl < 0.0 or d <= 0.0:
                continue
            diffusion_creep['V'] = Vdiff
            diffusion_creep['d'] = d
            dislocation_creep['V'] = Vdisl
            dislocation_creep['d'] = d

            # solve for A
            depth1 = 250e3
            T1 = T_func(depth1)
            P1 = P_func(depth1)
            # diffusion creep
            eta_diff1 = 2e20
            diff_A = CreepComputeA(diffusion_creep, strain_rate, P1, T1, eta_diff1)
            diffusion_creep['A'] = diff_A
            eta_disl1 = 2e20
            disl_A = CreepComputeA(dislocation_creep, strain_rate, P1, T1, eta_disl1, use_effective_strain_rate=True)
            dislocation_creep['A'] = disl_A

            # 660 km 
            # diffusion creep
            T660 = T_func(660e3)
            P660 = P_func(660e3)
            eta_diff660 = CreepRheology(diffusion_creep, strain_rate, P660, T660)
            # dislocation creep
            eta_disl660 = CreepRheology(dislocation_creep, strain_rate, P660, T660, use_effective_strain_rate=True)
            eta660 = ComputeComposite(eta_diff660, eta_disl660)
            
            # other constraints
            depth2 = 300e3
            Ttemp = T_func(depth2)
            Ptemp = P_func(depth2)
            # diffusion creep
            eta_diff2 = CreepRheology(diffusion_creep, strain_rate, Ptemp, Ttemp)
            # dislocation creep
            eta_disl2 = CreepRheology(dislocation_creep, strain_rate, Ptemp, Ttemp, use_effective_strain_rate=True)

            if eta660 >= eta660range[0] and eta660 <= eta660range[1] and eta_disl2 < eta_diff2:
                new_pair = {'diff': diffusion_creep, 'disl': dislocation_creep}
                constrained_rheologies.append(new_pair)
                constrained_ds.append(d)
                constrained_Vdisls.append(Vdisl)
                constrained_Vdiffs.append(Vdiff)
        
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(constrained_Vdiffs, constrained_Vdisls, constrained_ds)
        ax.set_xlabel('Vdiff [m^3/mol]')
        ax.set_ylabel('Vdisl [m^3/mol]')
        # ax.set_zlabel('d [um]')
        fig_name = 'constrained_rheology_%s_N%d.png' % (rheology, N)
        fig_path = os.path.join(RESULT_DIR, fig_name)
        fig.savefig(fig_path)


def Config(_kwargs, _name, _default):
    """
    def Config(_kwargs, _name, _default)

    Read variable value and assign default if not found
    """
    try:
        value = _kwargs[_name]
    except KeyError:
        value = _default
    return value


def CreepStress(creep_type, strain_rate, P, T, d, Coh):
    """
    def DislocationCreep(strain_rate, P, T, d, Coh)

    Calculate stress by flow law in form of (strain_rate / B)^(1.0 / n) * exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Mpa
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    # calculate B
    B = A * d**(-p) * Coh**r
    return (strain_rate / B)**(1.0 / n) * np.exp((E + P * V) / (n * R * T))


def CreepRheology(creep_type, strain_rate, P, T, d=None, Coh=None, **kwargs):
    """
    def CreepRheology(creep_type, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep_type['d']
    if Coh is None:
        Coh = creep_type['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    eta = F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6
    return eta


def CreepComputeA(creep_type, strain_rate, P, T, eta, d=None, Coh=None, **kwargs):
    """
    def CreepRheology(creep_type, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep_type['d']
    if Coh is None:
        Coh = creep_type['Coh']
    # calculate B
    B = (F/eta)**n * strain_rate**(1-n) * np.exp((E+P*V)/(R*T)) * (1e6)**n
    A = B * d**p * Coh**(-r)
    return A


def CreepRheologyInAspectViscoPlastic(creep_type, strain_rate, P, T):
    """
    def CreepRheologyInAspectVisoPlastic(creep_type, strain_rate, P, T)

    Calculate viscosity by way of Visco Plastic module in aspect
    flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: m
     - Return value: Pa*s
    """
    A = creep_type['A']
    m = creep_type['m']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    d = creep_type['d']
    # calculate B
    return 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))


def Convert2AspectInput(creep_type, **kwargs):
    """
    Viscosity is calculated by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6
    while in aspect, flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    Original Units:
     - P: Pa
     - T: K
     - d: mm
     - Coh: H / 10^6 Si
    Original Units:
     - P: Pa
     - T: K
     - d: m
    """
    # read in initial value
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    d = creep_type['d']
    Coh = creep_type['Coh']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0
    # prepare values for aspect
    aspect_creep_type = {}
    # stress in the original equation is in Mpa, grain size is in um
    aspect_creep_type['A'] = 1e6**(-p) * (2e6)**(-n) * Coh**r * A / F**n  # F term: use effective strain rate
    aspect_creep_type['d'] = d / 1e6
    aspect_creep_type['n'] = n
    aspect_creep_type['m'] = p
    aspect_creep_type['E'] = E
    aspect_creep_type['V'] = V
    return aspect_creep_type


def GetLowerMantleRheology(upper_mantle_creep_method, jump, T, P, **kwargs):
    """
    get flow law parameters in lower mantle based on upper mantle viscosity and jump in viscosity
    variables:
     - jump: viscosity jump at 660km
     - T: temperature at 660km
     - P: pressure at 660km
    """
    # extra inputs
    strategy = Config(kwargs, 'strategy', 'A')
    V1 = Config(kwargs, 'V1', upper_mantle_creep_method['V'])
    # read upper mantle values
    A = upper_mantle_creep_method['A']
    m = upper_mantle_creep_method['m']
    n = upper_mantle_creep_method['n']
    E = upper_mantle_creep_method['E']
    V = upper_mantle_creep_method['V']
    d = upper_mantle_creep_method['d']

    lower_mantle_creep_method = dict(upper_mantle_creep_method)
    lower_mantle_creep_method['V'] = V1
    if strategy is 'A':
        lower_mantle_creep_method['A'] = jump**(-n) * A * math.exp(P * (V1 - V) / (R * T))
    else:
        lower_mantle_creep_method['d'] = jump**(n / m) * d * math.exp(P * (V-V1) / (m * R * T))
    return lower_mantle_creep_method


def ComputeComposite(eta_diff, eta_disl):
    '''
    compute value of composite viscosity from value of diffusion creep and 
    dislocation creep.
    '''
    eta = 1.0 / (1.0/eta_diff + 1.0/eta_disl)
    return eta


def ReadAspectProfile(depth_average_path):
    """
    read a T,P profile from aspect's depth average file
    """
    # check file exist
    assert(os.access(depth_average_path, os.R_OK))
    # read that
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.ReadHeader(depth_average_path)
    DepthAverage.ReadData(depth_average_path)
    DepthAverage.SplitTimeStep()
    i0 = DepthAverage.time_step_indexes[0][-1] * DepthAverage.time_step_length
    i1 = DepthAverage.time_step_indexes[1][0] * DepthAverage.time_step_length
    data = DepthAverage.data[i0:i1, :]
    col_depth = DepthAverage.header['depth']['col']
    col_P = DepthAverage.header['adiabatic_pressure']['col']
    col_T = DepthAverage.header['temperature']['col']
    depths = data[:, col_depth]
    pressures = data[:, col_P]
    temperatures = data[:, col_T]
    return depths, pressures, temperatures


def PlotAlongProfile(depths, pressures, temperatures, fig_path_base, **kwargs):
    '''
    plot along a T, P profile in aspect
    '''
    # compute viscosity
    eta_annotation = '' # use this to annotate figure title
    rheology = kwargs.get('rheology', 'HK03')
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
    dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
    strain_rate = 1e-15

    # grain size
    d = kwargs.get('d', 1e4)
    
    # diffusion creep
    implementation = kwargs.get('implementation', 'LHY')
    if implementation == 'LHY':
        diffusion_creep['d'] = d
        eta_diff = CreepRheology(diffusion_creep, strain_rate, pressures, temperatures)
        eta_annotation += rheology
    elif implementation == 'MB':
        # use magali's implementation
        coh = 1000
        water = 'con' # 'wet, dry, con'
        mod = 'orig'
        Edev = 'mid'
        Vdev = 'mid'
        eta_diff = visc_diff_HK(temperatures,pressures,d,coh,water,mod,Edev,Vdev)
        eta_annotation += '%s_%s_E%s_V%s' % (water,mod, Edev, Vdev)
    else:
        raise CheckValueError('%s is not a valid implementation' % implementation)
    # dislocation creep
    dislocation_creep['d'] = d
    eta_disl = CreepRheology(dislocation_creep, strain_rate, pressures, temperatures, use_effective_strain_rate=True)
    eta = ComputeComposite(eta_diff, eta_disl)

    # plot
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    color = 'tab:blue'
    axs[0].plot(pressures/1e9, depths/1e3, color=color, label='pressure')
    axs[0].set_ylabel('Depth [km]') 
    axs[0].set_xlabel('Pressure [GPa]', color=color) 
    # axs[0].invert_yaxis()
    ylim=[660.0, 0.0]
    axs[0].set_ylim(ylim)
    # ax2: temperature
    color = 'tab:red'
    ax2 = axs[0].twiny()
    ax2.set_ylim(ylim)
    ax2.plot(temperatures, depths/1e3, color=color, label='temperature')
    ax2.set_xlabel('Temperature [K]', color=color) 
    # second: viscosity
    axs[1].semilogx(eta_diff, depths/1e3, 'c', label='diffusion creep')
    axs[1].semilogx(eta_disl, depths/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
    axs[1].semilogx(eta, depths/1e3, 'r--', label='Composite')
    axs[1].set_xlim([1e18,1e25])
    axs[1].set_ylim(ylim)
    # axs[1].invert_yaxis()
    axs[1].grid()
    axs[1].set_ylabel('Depth [km]') 
    axs[1].set_xlabel('Viscosity [Pa*s]')
    _title = 'Viscosity (%s)' % eta_annotation
    axs[1].set_title(_title)
    axs[1].legend()
    fig.tight_layout()
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s_%s_d%.2e_%s.%s" % (fig_path_base0, rheology, d, implementation, fig_path_type)
    plt.savefig(fig_path)
    print("New figure: %s" % fig_path)


def PlotAlongProfileJson(depths, pressures, temperatures, file_path, fig_path_base):
    '''
    plot along a T, P profile in aspect
    '''
    # compute viscosity
    eta_annotation = '' # use this to annotate figure title
    rheology = 'Aspect'
    eta_annotation += rheology
    strain_rate = 1e-15

    if (file_path.rpartition('.')[-1] == 'json'):
        with open(file_path, 'r') as fin:
            Rheology = json.load(fin)
        diffusion_creep = Rheology['diffusion_creep']
        dislocation_creep = Rheology['dislocation_creep']
        eta_annotation += '-json'
    elif (file_path.rpartition('.')[-1] == 'prm'):
        with open(file_path, 'r') as fin:
            inputs = ParseFromDealiiInput(fin)
        diffusion_creep, dislocation_creep = UpperMantleRheologyViscoPlastic(inputs)
        eta_annotation += '-prm'
    else:
        raise FileNotFoundError('Configuration file must be json or prm')

    # screen output 
    print('read rheology parameterization(diff, disl):')
    print(diffusion_creep)
    print(dislocation_creep)
    

    # diffusion creep
    eta_diff = CreepRheologyInAspectViscoPlastic(diffusion_creep, strain_rate, pressures, temperatures)
   
    # dislocation creep
    eta_disl = CreepRheologyInAspectViscoPlastic(dislocation_creep, strain_rate, pressures, temperatures)
    
    eta = ComputeComposite(eta_diff, eta_disl)

    # plot
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    color = 'tab:blue'
    axs[0].plot(pressures/1e9, depths/1e3, color=color, label='pressure')
    axs[0].set_ylabel('Depth [km]') 
    axs[0].set_xlabel('Pressure [GPa]', color=color) 
    # axs[0].invert_yaxis()
    ylim=[660.0, 0.0]
    axs[0].set_ylim(ylim)
    # ax2: temperature
    color = 'tab:red'
    ax2 = axs[0].twiny()
    ax2.set_ylim(ylim)
    ax2.plot(temperatures, depths/1e3, color=color, label='temperature')
    ax2.set_xlabel('Temperature [K]', color=color) 
    # second: viscosity
    axs[1].semilogx(eta_diff, depths/1e3, 'c', label='diffusion creep')
    axs[1].semilogx(eta_disl, depths/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
    axs[1].semilogx(eta, depths/1e3, 'r--', label='Composite')
    axs[1].set_xlim([1e18,1e25])
    axs[1].set_ylim(ylim)
    # axs[1].invert_yaxis()
    axs[1].grid()
    axs[1].set_ylabel('Depth [km]') 
    axs[1].set_xlabel('Viscosity [Pa*s]')
    _title = 'Viscosity (%s)' % eta_annotation
    axs[1].set_title(_title)
    axs[1].legend()
    fig.tight_layout()
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s_%s_d%.2e.%s" % (fig_path_base0, rheology, diffusion_creep['d'], fig_path_type)
    plt.savefig(fig_path)
    print("New figure: %s" % fig_path)


def ConstrainASPECT(file_path):
    '''
    Figure out contrain of rheology by random walk
    Inputs:
        file_path(str): a profile from ASPECT
    '''
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(file_path)
    # do a random walk
    Operator.ConstrainRheology()
    

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
    parser.add_argument('-j', '--json', type=str,
                        default=None,
                        help='path to a json file')
    parser.add_argument('-r', '--rheology', type=str,
                        default='HK03',
                        help='Type of rheology to use')
    parser.add_argument('-d', '--grain_size', type=float,
                        default=1e4,
                        help='Grain Size')
    parser.add_argument('-P', '--pressure', type=float,
                        default=10e9,
                        help='Pressure (Pa)')
    parser.add_argument('-T', '--temperature', type=float,
                        default=1673,
                        help='Temperature (K)')
    parser.add_argument('-S', '--strain_rate', type=float,
                        default=1e-15,
                        help='Strain Rate (s^-1)')
    parser.add_argument('-E', '--use_effective_strain_rate', type=int,
                        default=0,
                        help='If use effective strain rate instead of experimental value (0 or 1)')
    parser.add_argument('-im', '--implementation', type=str,
                        default='LHY',
                        help='implementation of rheology(LHY or MB)')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'convert_to_ASPECT':
        rheology = arg.rheology
        # read in standard flow law parameters
        RheologyPrm = RHEOLOGY_PRM()
        diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
        dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
        # convert 2 aspect
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=True)
        # save to output
        if arg.json is not None:
            creep_in_aspect = {}
            creep_in_aspect['diffusion_creep'] = diffusion_creep_aspect
            creep_in_aspect['dislocation_creep'] = dislocation_creep_aspect
            with open(arg.json, 'w') as fout:
                json.dump(creep_in_aspect, fout)
        # screen output
        print("ASPECT diffusion creep: ")
        print(diffusion_creep_aspect)
        print("ASPECT dislocation creep: ")
        print(dislocation_creep_aspect)
    
    elif _commend == 'compute_creep_viscosity':
        rheology = arg.rheology
        # read in standard flow law parameters
        RheologyPrm = RHEOLOGY_PRM()
        diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
        dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
        eta_diff = CreepRheology(diffusion_creep, arg.strain_rate, arg.pressure, arg.temperature)
        eta_disl = CreepRheology(dislocation_creep, arg.strain_rate, arg.pressure, arg.temperature, use_effective_strain_rate=arg.use_effective_strain_rate)
        # screen output
        print("eta_diff = %4e" % eta_diff)
        print("eta_disl = %4e" % eta_disl)

    elif _commend == 'compute_ASPECT_viscosity':
        # read from json file
        with open(arg.json, 'r') as fin:
            Rheology = json.load(fin)
        diffusion_creep = Rheology['diffusion_creep']
        dislocation_creep = Rheology['dislocation_creep']
        eta_diff = CreepRheologyInAspectViscoPlastic(diffusion_creep, arg.strain_rate, arg.pressure, arg.temperature)
        eta_disl = CreepRheologyInAspectViscoPlastic(dislocation_creep, arg.strain_rate, arg.pressure, arg.temperature)
        # screen output
        print("eta_diff = %4e" % eta_diff)
        print("eta_disl = %4e" % eta_disl)
    
    elif _commend == 'plot_along_aspect_profile':
        # todo
        fig_path = os.path.join(RESULT_DIR, 'along_profile_rhoelogy.png')
        depths, pressures, temperatures = ReadAspectProfile(arg.inputs)
        PlotAlongProfile(depths, pressures, temperatures, fig_path, rheology=arg.rheology, d=arg.grain_size, implementation=arg.implementation)
    
    elif _commend == 'plot_along_aspect_profile_with_json':
        fig_path = os.path.join(RESULT_DIR, 'along_profile_rhoelogy.png')
        depths, pressures, temperatures = ReadAspectProfile(arg.inputs)
        PlotAlongProfileJson(depths, pressures, temperatures, arg.json, fig_path)

    elif _commend == 'constrain_aspect_rheology':
        ConstrainASPECT(arg.inputs)
    
    else:
        raise CheckValueError('%s is not a valid commend' % _commend)



# run script
if __name__ == '__main__':
    main()
