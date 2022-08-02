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

        python -m shilofue.Rheology plot_along_aspect_profile -r HK03 -i files/TwoDSubduction/reference/depth_average.txt -im MB

        -im : LHY(default) or MB, use my implementation of formula or Magali's (in script flow_law_function)
  
  - plot along a profile from aspect, using rheology parameterization from a prm file to check implementation in code

        python -m shilofue.Rheology plot_along_aspect_profile_with_json 
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/output/depth_average.txt 
        -j /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/case.prm

        -i: depth_average file
        -j: prm file / json file, in case of a json file, it contains diffusion + dislocation creep(as output from the convert_to_ASPECT command)

  - constrain rheology by random sampling the parameter space.

        python -m shilofue.Rheology constrain_aspect_rheology 
        -i files/TwoDSubduction/reference/depth_average.txt 
        -v 1
        --save_profile 1 --include_lower_mantle 30.0

        -v: version, default is 0
        --save_profile: save separate accepted profile to png files
        --include_lower_mantle: jump between lower mantle / upper mantle rheology

  - hand modify the rheology parameters and see the plotting

        python -m shilofue.Rheology derive_mantle_rheology -i files/TwoDSubduction/reference/depth_average.txt  -v 1 --save_profile 1 -r HK03_wet_mod --diff "0.33333333333,-40e3,-5.5e-6,1.73205080757,20e3,0.0"

descriptions
""" 

import json
import os
import sys
import math
import argparse
import numpy as np
from scipy.special import erf
from matplotlib import pyplot as plt
from matplotlib import gridspec
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
from shilofue.FlowLaws import visc_diff_HK
from shilofue.ParsePrm import ParseFromDealiiInput, UpperMantleRheologyViscoPlastic
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from shutil import rmtree

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
        # this is found in the supplementary material in the paper.
        # Note that the original value in the paper uses "Pa" for A,
        # while I converted it to "MPa" here.
        self.AB17_disl = \
            {
                "A": 25.7,
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
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        # 'wet' indicates this has to applied with a rheology of water
        self.HK03_wet_mod_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3,
                "V" : 23e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.5,
                "E" : 520e3,
                "V" : 24e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # I bring the values to the limit of the range
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        self.HK03_wet_mod1_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3 - 25e3,
                "V" : 23e-6 -5.5e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod1_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.5,
                "E" : 520e3 + 40e3,
                "V" : 24e-6 + 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # I bring the values to the limit of the range
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        # this is specifically the one I used for the TwoD models.
        # Combined with the usage of function "MantleRheology_v0", then same rheology
        # could be reproduced
        self.HK03_wet_mod_2d_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3 - 40e3,
                "V" : 23e-6 -5.5e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod_2d_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.5,
                "E" : 520e3 + 20e3,
                "V" : 24e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        
        self.water = \
            {
                "A" : 87.75,             # H/(10^6*Si)/MPa
                "E" : 50e3,                     # J/mol +/-2e3
                "V" : 10.6e-6                     # m^3/mol+/-1
            }

        # this is the values used in the ARCAY17 paper
        # note: their rheology is only stress dependent (dislocation creep)
        # their yielding criterion is stress dependent as well.
        self.ARCAY17_diff = None
        self.ARCAY17_disl = \
            {
                "A" : 339428.7,
                "p" : 0.0,
                "r" : 0.0,  # not dependent on the "Coh"
                "n" : 3.0,
                "E" : 465e3,
                "V" : 17e-6,
                "d" : 1e4, # not dependent on d
                "Coh" : 1000.0
            }
        self.ARCAY17_plastic = \
        {
            "friction" : 0.05,
            "cohesion": 1e6, # pa
            "n": 30.0,
            "ref strain rate" : 1.0e-14,
            "type": "stress dependent"
        }

        
        
        self.water = \
            {
                "A" : 87.75,             # H/(10^6*Si)/MPa
                "E" : 50e3,                     # J/mol +/-2e3
                "V" : 10.6e-6                     # m^3/mol+/-1
            }

    def get_rheology(self, _name, _type):
        '''
        read rheology parameters, and account for effects of water if it is a wet rheology
        '''
        assert(_type in ['diff', 'disl', 'plastic'])
        _attr = _name + "_" + _type
        if not hasattr(self, _attr):
            raise ValueError("RHEOLOGY_PRM object doesn't have attribute %s" % _attr)
        creep = getattr(self, _attr)
        if "wet" in creep:
            # foh enters explicitly, converting to use Coh
            assert(_type in ['diffusion', 'dislocation'])
            ### effects of water accounted, see Magali's file explain_update_modHK03_rheology eq(5)
            water_creep = getattr(RheologyPrm, "water")
            creep['A'] = creep['A'] / (water_creep['A'] ** creep['r'])
            creep['V'] = creep['V'] - water_creep['V'] * creep['r']
            creep['E'] = creep['E'] - water_creep['E'] * creep['r']
        return creep


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
        self.output_profile = None # for the figure plotted
        self.output_json = None
        self.output_aspect_json = None
        self.diff = None
        self.disl = None
        self.diff_type = None
        self.disl_type = None
        self.plastic = None
        pass

    def SetRheology(self, **kwargs):
        '''
        set rheology type with instances of rheology (i.e. a dictionary)
        '''
        self.diff = kwargs.get('diff', None)
        self.disl = kwargs.get('disl', None)
        self.plastic = kwargs.get('plastic', None)
        pass
    
    def SetRheologyByName(self, **kwargs):
        '''
        set rheology type with instances of rheology (i.e. a dictionary)
        '''
        diff_type = kwargs.get('diff', None)
        disl_type = kwargs.get('disl', None)
        self.diff_type = diff_type
        self.disl_type = disl_type
        plastic_type = kwargs.get('plastic', None)
        if diff_type != None:
            self.diff = self.RheologyPrm.get_rheology(diff_type, 'diff')
        if disl_type != None:
            self.disl = self.RheologyPrm.get_rheology(disl_type, 'disl')
        if plastic_type != None:
            self.plastic = self.RheologyPrm.get_rheology(plastic_type, 'plastic')

    
    def ReadProfile(self, file_path):
        self.depths, self.pressures, self.temperatures = ReadAspectProfile(file_path)

    def MantleRheology_v0(self, **kwargs):
        '''
        Derive mantle rheology from an aspect profile
        '''
        use_effective_strain_rate = kwargs.get('use_effective_strain_rate', True)
        strain_rate = kwargs.get('strain_rate', 1e-15)
        eta_diff = np.ones(self.depths.size)
        eta_disl = np.ones(self.depths.size)
        eta = np.ones(self.depths.size)
        dEdiff = kwargs.get('dEdiff', 0.0)  # numbers for the variation in the rheology
        dVdiff = kwargs.get('dVdiff', 0.0)
        dEdisl = kwargs.get('dEdisl', 0.0)
        dVdisl = kwargs.get('dVdisl', 0.0)
        rheology = kwargs.get('rheology', 'HK03_wet_mod')
        save_profile = kwargs.get('save_profile', 0)
        save_json = kwargs.get('save_json', 0)
        
        
        diffusion_creep, dislocation_creep = GetRheology(rheology)
        diffusion_creep['E'] += dEdiff
        dislocation_creep['E'] += dEdisl
        diffusion_creep['V'] += dVdiff
        dislocation_creep['V'] += dVdisl

        # convert T, P as function
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # < 410 km
        depth_up = 410e3
        depth_low = 660e3
        mask_up = (self.depths <= depth_up)
        eta_diff[mask_up] = CreepRheology(diffusion_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up])
        eta_disl[mask_up] = CreepRheology(dislocation_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=use_effective_strain_rate)
        eta[mask_up] = ComputeComposite(eta_diff[mask_up], eta_disl[mask_up])


        # MTZ
        mask_mtz = (self.depths > depth_up) & (self.depths < depth_low)
        if True:
            # MTZ from olivine rheology
            eta_diff[mask_mtz] = CreepRheology(diffusion_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz])
            eta_disl[mask_mtz] = CreepRheology(dislocation_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=use_effective_strain_rate)
            eta[mask_mtz] = ComputeComposite(eta_diff[mask_mtz], eta_disl[mask_mtz])
        
        # lower mantle
        # diffusion creep
        mask_low = (self.depths > depth_low)
        jump_lower_mantle = kwargs.get('jump_lower_mantle', 30.0)
        # values for computing V
        depth_lm = 660e3
        depth_for_lm_mean = min(np.max(self.depths), 1700e3)
        T_lm_mean = T_func(depth_for_lm_mean)
        P_lm_mean = P_func(depth_for_lm_mean)
        depth_max = self.depths[-1] - 10e3
        lm_grad_T = (T_func(depth_max) - T_func(depth_lm)) / (depth_max - depth_lm)
        lm_grad_P = (P_func(depth_max) - P_func(depth_lm)) / (depth_max - depth_lm)
        T660 = T_func(depth_lm)
        P660 = P_func(depth_lm)
        eta_diff660 = CreepRheology(diffusion_creep, strain_rate, P660, T660)
        # dislocation creep
        eta_disl660 = CreepRheology(dislocation_creep, strain_rate, P660, T660, use_effective_strain_rate=use_effective_strain_rate)
        eta660 = ComputeComposite(eta_diff660, eta_disl660)
        diff_lm = diffusion_creep.copy()
        # diff_lm['V'] = LowerMantleV(diffusion_creep['E'], T_lm_mean, P_lm_mean, lm_grad_T, lm_grad_P) # compute from less variation criteria
        diff_lm['V'] = 3e-6  # assign a value
        diff_lm['A'] = CreepComputeA(diff_lm, strain_rate, P660, T660, eta660*jump_lower_mantle)
        eta_diff[mask_low] = CreepRheology(diff_lm, strain_rate, self.pressures[mask_low], self.temperatures[mask_low])
        eta_disl[mask_low] = None  # this is just for visualization
        eta[mask_low] = eta_diff[mask_low]  # diffusion creep is activated in lower mantle

        # haskel constraint
        radius = 6371e3
        lith_depth = 100e3
        integral_depth = 1400e3
        mask_integral = (self.depths > lith_depth) & (self.depths < integral_depth)
        integral_cores = 4 * np.pi * (radius - self.depths)**2.0
        # upper mantle
        # use harmonic average
        # lower mantle
        integral = np.trapz(integral_cores[mask_integral] * np.log10(eta[mask_integral]), self.depths[mask_integral])
        volume = np.trapz(integral_cores[mask_integral], self.depths[mask_integral])
        average_log_eta = integral / volume
        # dump json file 
        constrained_rheology = {'diffusion_creep': diffusion_creep, 'dislocation_creep': dislocation_creep, 'diffusion_lm': diff_lm}
        # convert aspect rheology
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
        diffusion_lm_aspect = Convert2AspectInput(diff_lm)
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=use_effective_strain_rate)
        constrained_rheology_aspect = {'diffusion_creep': diffusion_creep_aspect, 'dislocation_creep': dislocation_creep_aspect, 'diffusion_lm': diffusion_lm_aspect}
        if save_json == 1:
            json_path = os.path.join(RESULT_DIR, "mantle_profile_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            json_path_aspect = os.path.join(RESULT_DIR, "mantle_profile_aspect_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            with open(json_path, 'w') as fout:
                json.dump(constrained_rheology, fout)
            with open(json_path_aspect, 'w') as fout:
                json.dump(constrained_rheology_aspect, fout)
            print("New json: %s" % json_path)
            print("New json: %s" % json_path_aspect)
            self.output_json = json_path
            self.output_json_aspect = json_path_aspect
        # plot
        if save_profile == 1:
            # plots
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            color = 'tab:blue'
            axs[0].plot(self.pressures/1e9, self.depths/1e3, color=color, label='pressure')
            axs[0].set_ylabel('Depth [km]') 
            axs[0].set_xlabel('Pressure [GPa] P660: %.4e' % (P660), color=color) 
            # axs[0].invert_yaxis()
            ylim=[np.ceil(np.max(self.depths) / 100e3) * 100.0, 0.0]
            axs[0].set_ylim(ylim)
            # ax2: temperature
            color = 'tab:red'
            ax2 = axs[0].twiny()
            ax2.set_ylim(ylim)
            ax2.plot(self.temperatures, self.depths/1e3, color=color, label='temperature')
            ax2.set_xlabel('Temperature [K] T660: %.4e' % (T660), color=color) 
            # second: viscosity
            #   upper mantle
            axs[1].semilogx(eta_diff, self.depths/1e3, 'c', label='diffusion creep')
            axs[1].semilogx(eta_disl, self.depths/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
            axs[1].semilogx(eta, self.depths/1e3, 'r--', label='Composite')
            axs[1].set_xlim([1e18,1e25])
            axs[1].set_ylim(ylim)
            # axs[1].invert_yaxis()
            axs[1].grid()
            axs[1].set_ylabel('Depth [km]')
            axs[1].legend()
            axs[1].set_title('%s_lowerV_%.4e_haskell%.2f' % (rheology, diff_lm['V'], average_log_eta))
            # save figure
            fig_path = os.path.join(RESULT_DIR, "mantle_profile_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.png" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            fig.savefig(fig_path)
            print("New figure: %s" % fig_path)
            plt.close()
            self.output_profile = fig_path
            pass
        return constrained_rheology_aspect
    
    def MantleRheology_v1(self, **kwargs):
        '''
        Derive mantle rheology from an aspect profile
        In this version, I would use the F factor as the default for computing the viscosity.
        Also I am going to use CreepRheology_v1_v1 function rather thant the CreepRheology_v1 function,
        where the factor of F is dealt with correctly.
        '''
        strain_rate = kwargs.get('strain_rate', 1e-15)
        eta_diff = np.ones(self.depths.size)
        eta_disl = np.ones(self.depths.size)
        eta_disl13 = np.ones(self.depths.size)
        eta13 = np.ones(self.depths.size)
        eta = np.ones(self.depths.size)
        dEdiff = float(kwargs.get('dEdiff', 0.0))  # numbers for the variation in the rheology
        dVdiff = float(kwargs.get('dVdiff', 0.0))
        dAdiff_ratio = float(kwargs.get("dAdiff_ratio", 1.0))
        dAdisl_ratio = float(kwargs.get("dAdisl_ratio", 1.0))
        dEdisl = float(kwargs.get('dEdisl', 0.0))
        dVdisl = float(kwargs.get('dVdisl', 0.0))
        rheology = kwargs.get('rheology', 'HK03_wet_mod')
        save_profile = kwargs.get('save_profile', 0)
        save_json = kwargs.get('save_json', 0)
        
        
        diffusion_creep, dislocation_creep = GetRheology(rheology)
        diffusion_creep['A'] *= dAdiff_ratio
        diffusion_creep['E'] += dEdiff
        dislocation_creep['A'] *= dAdisl_ratio
        dislocation_creep['E'] += dEdisl
        diffusion_creep['V'] += dVdiff
        dislocation_creep['V'] += dVdisl

        # convert T, P as function
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # < 410 km
        depth_up = 410e3
        depth_low = 660e3
        mask_up = (self.depths < depth_up)
        eta_diff[mask_up] = CreepRheology_v1(diffusion_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=True)
        eta_disl[mask_up] = CreepRheology_v1(dislocation_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=True)
        eta_disl13[mask_up] = CreepRheology_v1(dislocation_creep, 1e-13, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=True)
        eta[mask_up] = ComputeComposite(eta_diff[mask_up], eta_disl[mask_up])
        eta13[mask_up] = ComputeComposite(eta_diff[mask_up], eta_disl13[mask_up])


        # MTZ
        mask_mtz = (self.depths > depth_up) & (self.depths < depth_low)
        if True:
            # MTZ from olivine rheology
            eta_diff[mask_mtz] = CreepRheology_v1(diffusion_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=True)
            eta_disl[mask_mtz] = CreepRheology_v1(dislocation_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=True)
            eta_disl13[mask_mtz] = CreepRheology_v1(dislocation_creep, 1e-13, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=True)
            eta[mask_mtz] = ComputeComposite(eta_diff[mask_mtz], eta_disl[mask_mtz])
            eta13[mask_mtz] = ComputeComposite(eta_diff[mask_mtz], eta_disl13[mask_mtz])
        
        # lower mantle
        # diffusion creep
        mask_low = (self.depths > depth_low)
        jump_lower_mantle = kwargs.get('jump_lower_mantle', 30.0)
        # values for computing V
        depth_lm = 660e3
        T_lm_mean = T_func(1700e3)
        P_lm_mean = P_func(1700e3)
        depth_max = self.depths[-1] - 10e3
        lm_grad_T = (T_func(depth_max) - T_func(depth_lm)) / (depth_max - depth_lm)
        lm_grad_P = (P_func(depth_max) - P_func(depth_lm)) / (depth_max - depth_lm)
        T660 = T_func(depth_lm)
        P660 = P_func(depth_lm)
        eta_diff660 = CreepRheology_v1(diffusion_creep, strain_rate, P660, T660, use_effective_strain_rate=True)
        # dislocation creep
        eta_disl660 = CreepRheology_v1(dislocation_creep, strain_rate, P660, T660, use_effective_strain_rate=True)
        eta660 = ComputeComposite(eta_diff660, eta_disl660)
        diff_lm = diffusion_creep.copy()
        # diff_lm['V'] = LowerMantleV(diffusion_creep['E'], T_lm_mean, P_lm_mean, lm_grad_T, lm_grad_P) # compute from less variation criteria
        diff_lm['V'] = 3e-6  # assign a value
        diff_lm['A'] = CreepComputeA_v1(diff_lm, strain_rate, P660, T660, eta660*jump_lower_mantle, use_effective_strain_rate=True)
        # print('diff_lm: ', diff_lm) # print additional information
        eta_diff[mask_low] = CreepRheology_v1(diff_lm, strain_rate, self.pressures[mask_low], self.temperatures[mask_low], use_effective_strain_rate=True)
        eta_disl[mask_low] = None  # this is just for visualization
        eta_disl13[mask_low] = None  # this is just for visualization
        eta[mask_low] = eta_diff[mask_low]  # diffusion creep is activated in lower mantle
        eta13[mask_low] = eta_diff[mask_low]  # diffusion creep is activated in lower mantle

        # haskel constraint
        radius = 6371e3
        lith_depth = 100e3
        integral_depth = 1400e3
        mask_integral = (self.depths > lith_depth) & (self.depths < integral_depth)
        integral_cores = 4 * np.pi * (radius - self.depths)**2.0
        # upper mantle
        # use harmonic average
        # lower mantle
        integral = np.trapz(integral_cores[mask_integral] * np.log10(eta[mask_integral]), self.depths[mask_integral])
        volume = np.trapz(integral_cores[mask_integral], self.depths[mask_integral])
        average_log_eta = integral / volume
        # dump json file 
        constrained_rheology = {'diffusion_creep': diffusion_creep, 'dislocation_creep': dislocation_creep, 'diffusion_lm': diff_lm}
        # convert aspect rheology
        diffusion_creep_aspect = Convert2AspectInput_v1(diffusion_creep, use_effective_strain_rate=True)
        diffusion_lm_aspect = Convert2AspectInput_v1(diff_lm, use_effective_strain_rate=True)
        dislocation_creep_aspect = Convert2AspectInput_v1(dislocation_creep, use_effective_strain_rate=True)
        constrained_rheology_aspect = {'diffusion_creep': diffusion_creep_aspect, 'dislocation_creep': dislocation_creep_aspect, 'diffusion_lm': diffusion_lm_aspect}
        if save_json == 1:
            json_path = os.path.join(RESULT_DIR, "mantle_profile_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            json_path_aspect = os.path.join(RESULT_DIR, "mantle_profile_aspect_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            with open(json_path, 'w') as fout:
                json.dump(constrained_rheology, fout)
            with open(json_path_aspect, 'w') as fout:
                json.dump(constrained_rheology_aspect, fout)
            print("New json: %s" % json_path)
            print("New json: %s" % json_path_aspect)
        # plot
        if save_profile == 1:
            # plots
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            color = 'tab:blue'
            axs[0].plot(self.pressures/1e9, self.depths/1e3, color=color, label='pressure')
            axs[0].set_ylabel('Depth [km]') 
            axs[0].set_xlabel('Pressure [GPa] P660: %.4e' % (P660), color=color) 
            # axs[0].invert_yaxis()
            ylim=[2890, 0.0]
            axs[0].set_ylim(ylim)
            # ax2: temperature
            color = 'tab:red'
            ax2 = axs[0].twiny()
            ax2.set_ylim(ylim)
            ax2.plot(self.temperatures, self.depths/1e3, color=color, label='temperature')
            ax2.set_xlabel('Temperature [K] T660: %.4e' % (T660), color=color) 
            # second: viscosity
            #   upper mantle
            axs[1].semilogx(eta_diff, self.depths/1e3, 'c', label='diffusion creep')
            axs[1].semilogx(eta_disl, self.depths/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
            axs[1].semilogx(eta, self.depths/1e3, 'r--', label='Composite')
            axs[1].set_xlim([1e18,1e25])
            axs[1].set_ylim(ylim)
            axs[1].grid()
            axs[1].set_ylabel('Depth [km]')
            axs[1].legend()
            axs[1].set_title('%s_lowerV_%.4e_haskell%.2f' % (rheology, diff_lm['V'], average_log_eta))
            # third, viscosity at 1e-13 /s strain rate
            axs[2].semilogx(eta_diff, self.depths/1e3, 'c', label='diffusion creep')
            axs[2].semilogx(eta_disl13, self.depths/1e3, 'g', label='dislocation creep(%.2e)' % 1e-13)
            axs[2].semilogx(eta13, self.depths/1e3, 'r--', label='Composite')
            axs[2].set_xlim([1e18,1e25])
            axs[2].set_ylim(ylim)
            axs[2].grid()
            axs[2].set_ylabel('Depth [km]')
            axs[2].legend()
            axs[2].set_title('strain_rate1.0e-13')
            # save figure
            fig_path = os.path.join(RESULT_DIR, "mantle_profile_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e.png" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl))
            fig.savefig(fig_path)
            print("New figure: %s" % fig_path)
            plt.close()
            pass
        return constrained_rheology_aspect
    
    def ConstrainRheology_v0(self, **kwargs):
        '''
        varying around a give rheology with variation with applied constraints
        Version 0:
            1. both values of viscosity at a depth of 250e3 equal 2e20(to compute A), 
            2. dislocation creep is smaller than diffusion creep at 300e3.  
            Lower mantle is then assigned so that: 
                1. only diffusion creep is activated. 
                2 viscosity is nearly constant by controlling V. 
                3. There is a 30 times jump between u/l boundary by controlling A. 
        inputs: 
            kwargs(dict):
                rheology: type of initial rheology
        '''
        constrained_rheologies = []
        constrained_ds = []
        constrained_Vdisls = []
        constrained_Vdiffs = []
        rheology = kwargs.get('rheology', 'HK03')
        # get rheology
        diffusion_creep, dislocation_creep = GetRheology('rheology')
        
        lm_diff = {}
        strain_rate = 1e-15
        radius = kwargs.get('radius', 6371e3)  # earth radius
        
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # lower mantle
        depth_lm = 660e3
        T_lm_mean = T_func(1700e3)
        P_lm_mean = P_func(1700e3)
        depth_max = self.depths[-1] - 10e3
        lm_grad_T = (T_func(depth_max) - T_func(depth_lm)) / (depth_max - depth_lm)
        lm_grad_P = (P_func(depth_max) - P_func(depth_lm)) / (depth_max - depth_lm)

        # grain size
        d_mean = kwargs.get('d', 0.75e4)

        Vdiff_sigma = 2e-6
        Vdisl_sigma = 5.5e-6
        d_sigma = 5e3

        # random walk
        N = 1001
        Vdiffs = np.random.normal(diffusion_creep['V'], Vdiff_sigma, N)
        Vdisls = np.random.normal(dislocation_creep['V'], Vdisl_sigma, N)
        ds = np.random.normal(d_mean, d_sigma, N)
        include_lower_mantle = kwargs.get('include_lower_mantle', None)
        mask_um = (self.depths < depth_lm)  # a mask to get the components of the upper mantle
        mask_lm = (self.depths >= depth_lm)  # a mask to get the components of the lower mantle
        eta660range = [4e20, 1e21]
        average_range = [0.65e21, 1.1e21]
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
            diff_A = CreepComputeA(diffusion_creep, strain_rate, P1, T1, eta_diff1) # key word 'd' is not used
            diffusion_creep['A'] = diff_A
            eta_disl1 = 2e20
            disl_A = CreepComputeA(dislocation_creep, strain_rate, P1, T1, eta_disl1, use_effective_strain_rate=True)
            dislocation_creep['A'] = disl_A

            # 660 km 
            # diffusion creep
            T660 = T_func(depth_lm)
            P660 = P_func(depth_lm)
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

            if include_lower_mantle is not None:
                # lower mantle rheology
                diff_lm = diffusion_creep.copy()
                diff_lm['V'] = LowerMantleV(diffusion_creep['E'], T_lm_mean, P_lm_mean, lm_grad_T, lm_grad_P)
                diff_lm['A'] = CreepComputeA(diff_lm, strain_rate, P660, T660, eta660*include_lower_mantle)

            
            # 1000km integral
            eta_diff = CreepRheology(diffusion_creep, strain_rate, self.pressures, self.temperatures)
            eta_disl = CreepRheology(dislocation_creep, strain_rate, self.pressures, self.temperatures, use_effective_strain_rate=True)
            eta = ComputeComposite(eta_diff, eta_disl)
            lith_depth = 100e3
            integral_depth = 1400e3
            mask_integral = (self.depths > lith_depth) & (self.depths < integral_depth)
            integral_cores = 4 * np.pi * (radius - self.depths)**2.0
            # upper mantle
            mask_um_integral = (mask_um & mask_integral)
            integral_um = np.trapz(eta[mask_um_integral] * integral_cores[mask_um_integral], self.depths[mask_um_integral])
            # lower mantle
            integral_lm = 0.0
            if include_lower_mantle is not None:
                mask_lm_integral = (mask_lm & mask_integral)
                eta_diff_lm = CreepRheology(diff_lm, strain_rate, self.pressures, self.temperatures)
                integral_lm = np.trapz(eta_diff_lm[mask_lm_integral] * integral_cores[mask_lm_integral], self.depths[mask_lm_integral])
            else:
                integral_lm = 4.0 / 3 * np.pi * ((radius - depth_lm)**3.0 - (radius - integral_depth)**3.0) * eta660 * 30.0 # assume 30 times jump
            volume = 4.0 / 3 * np.pi * ((radius - lith_depth)**3.0 - (radius - integral_depth)**3.0)
            average_eta = (integral_um + integral_lm) / volume

            conds = [(eta660 >= eta660range[0]) and (eta660 <= eta660range[1]), (eta_disl2 > eta_diff2),\
            (average_eta >= average_range[0]) and (average_eta <= average_range[1])]

            condition_indexes = [0, 1]
            cond_combined = True
            for i in condition_indexes:
                cond_combined = (cond_combined and conds[i])
            
            if eta660 >= eta660range[0] and eta660 <= eta660range[1] and eta_disl2 < eta_diff2:
                constrained_ds.append(d)
                constrained_Vdisls.append(Vdisl)
                constrained_Vdiffs.append(Vdiff)
                if include_lower_mantle is not None:
                    constrained_rheologies.append({'diff': diffusion_creep.copy(), 'disl': dislocation_creep.copy(), 'diff_lm': diff_lm.copy(),\
                    'average_upper_region': average_eta})
                else:
                    constrained_rheologies.append({'diff': diffusion_creep.copy(), 'disl': dislocation_creep.copy(),\
                    'average_upper_region': average_eta})

        # make a new directory
        fig_dir = os.path.join(RESULT_DIR, 'constrained_rheology_%s_N%d' % (rheology, N))
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        else:
            rmtree(fig_dir)
            os.mkdir(fig_dir)

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(constrained_Vdiffs, constrained_Vdisls, constrained_ds)
        ax.set_xlabel('Vdiff [m^3/mol]')
        ax.set_ylabel('Vdisl [m^3/mol]')
        ax.set_zlabel('d [um]')
        fig_name = 'constrained_rheology_%s_N%d.png' % (rheology, N)
        fig_path = os.path.join(fig_dir, fig_name)
        fig.savefig(fig_path)
        print("New figure: %s" % fig_path)
        plt.close()

        # plot profiles
        save_profile = kwargs.get('save_profile', 0)
        i = 0  # index
        for constrained_rheology in constrained_rheologies:
            # dump json
            Vdiff = constrained_rheology['diff']['V']
            Vdisl = constrained_rheology['disl']['V']
            d = constrained_rheology['disl']['d']
            json_path = os.path.join(fig_dir, "constrained_profile_Vdiff%.4e_Vdisl%.4e_d%.4e.json" % (Vdiff, Vdisl, d))
            with open(json_path, 'w') as fout:
                json.dump(constrained_rheology, fout)
                print("[%d / %d], New json: %s" % (i, len(constrained_rheologies), json_path))
            # save profile
            if save_profile == 1:
                eta_diff = CreepRheology(constrained_rheology['diff'], strain_rate, self.pressures, self.temperatures)
                eta_disl = CreepRheology(constrained_rheology['disl'], strain_rate, self.pressures, self.temperatures, use_effective_strain_rate=True)
                eta = ComputeComposite(eta_diff, eta_disl)
                # plots
                fig, axs = plt.subplots(1, 2, figsize=(10, 5))
                color = 'tab:blue'
                axs[0].plot(self.pressures/1e9, self.depths/1e3, color=color, label='pressure')
                axs[0].set_ylabel('Depth [km]') 
                axs[0].set_xlabel('Pressure [GPa]', color=color) 
                # axs[0].invert_yaxis()
                if include_lower_mantle is None:
                    ylim=[660.0, 0.0]
                else:
                    ylim=[2890, 0.0]
                axs[0].set_ylim(ylim)
                # ax2: temperature
                color = 'tab:red'
                ax2 = axs[0].twiny()
                ax2.set_ylim(ylim)
                ax2.plot(self.temperatures, self.depths/1e3, color=color, label='temperature')
                ax2.set_xlabel('Temperature [K]', color=color) 
                # second: viscosity
                #   upper mantle
                axs[1].semilogx(eta_diff[mask_um], self.depths[mask_um]/1e3, 'c', label='diffusion creep')
                axs[1].semilogx(eta_disl[mask_um], self.depths[mask_um]/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
                axs[1].semilogx(eta[mask_um], self.depths[mask_um]/1e3, 'r--', label='Composite')
                axs[1].set_xlim([1e18,1e25])
                axs[1].set_ylim(ylim)
                # axs[1].invert_yaxis()
                axs[1].grid()
                axs[1].set_ylabel('Depth [km]')
                axs[1].legend()
                if include_lower_mantle:
                    # include lower mantle info as title:
                    diff_lm = constrained_rheology['diff_lm']
                    eta_diff_lm = CreepRheology(constrained_rheology['diff_lm'], strain_rate, self.pressures, self.temperatures)
                    axs[1].semilogx(eta_diff_lm[mask_lm], self.depths[mask_lm]/1e3, 'c')
                    title_str = "lm_jump%.2d_Vlm%.4e_avg%.4e" % (include_lower_mantle, diff_lm['V'], constrained_rheology['average_upper_region'])
                    axs[1].set_title(title_str)
                # save figure
                fig_path = os.path.join(fig_dir, "constrained_profile_Vdiff%.4e_Vdisl%.4e_d%.4e.png" % (Vdiff, Vdisl, d))
                fig.savefig(fig_path)
                print("[%d / %d], New figure: %s" % (i, len(constrained_rheologies), fig_path))
                plt.close()
            i = i + 1
    
    def ConstrainRheology_v1(self, **kwargs):
        '''
        varying around a give rheology with variation with applied constraints
        Version 1:
            0. use modified wet rheology see the file from Magali
            1. compute V with rheology at 250km depth, this is in turn, constraint in a range
            2. dislocation creep is bigger than diffusion creep at 300e3.  
            3. value of rheology on 660 within a range
            Lower mantle is then assigned so that: 
                1. only diffusion creep is activated. 
                2 viscosity is nearly constant by controlling V. 
                3. There is a 30 times jump between u/l boundary by controlling A. 
        inputs: 
            kwargs(dict):
                rheology: type of initial rheology
        '''
        constrained_rheologies = []
        constrained_ds = []
        constrained_Vdisls = []
        constrained_Edisls = []
        constrained_Vdiffs = []
        constrained_Ediffs = []
        rheology = kwargs.get('rheology', 'HK03_wet_mod')
        N = 1001
        
        # make a new directory
        fig_dir = os.path.join(RESULT_DIR, 'constrained_rheology_v1_%s_N%d' % (rheology, N))
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        else:
            rmtree(fig_dir)
            os.mkdir(fig_dir)
        
        # get rheology
        diffusion_creep, dislocation_creep = GetRheology(rheology)
        diff_orig = diffusion_creep.copy()
        disl_orig = dislocation_creep.copy()
        orig_rheology = {'diff': diff_orig, 'disl': disl_orig}
        json_path = os.path.join(fig_dir, "original_profile.json")
        with open(json_path, 'w') as fout:
            json.dump(orig_rheology, fout)
            print("New json: %s" % json_path)
        
        lm_diff = {}
        strain_rate = 1e-15
        radius = kwargs.get('radius', 6371e3)  # earth radius
        
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # lower mantle
        depth_lm = 660e3
        T_lm_mean = T_func(1700e3)
        P_lm_mean = P_func(1700e3)
        depth_max = self.depths[-1] - 10e3
        lm_grad_T = (T_func(depth_max) - T_func(depth_lm)) / (depth_max - depth_lm)
        lm_grad_P = (P_func(depth_max) - P_func(depth_lm)) / (depth_max - depth_lm)

        # grain size
        d_mean = kwargs.get('d', 0.75e4)

        # range of sampling
        Vdiff_sigma = 5.5e-6
        Ediff_sigma = 40e3
        Vdisl_sigma = 4e-6 # 10e-6
        Edisl_sigma = 50e3 # 10e-6
        d_sigma = 5e3

        # random sampling
        Ediffs = np.random.normal(diffusion_creep['E'], Ediff_sigma, N)
        Edisls = np.random.normal(dislocation_creep['E'], Edisl_sigma, N)
        ds = np.random.normal(d_mean, d_sigma, N)
        include_lower_mantle = kwargs.get('include_lower_mantle', None)
        mask_um = (self.depths < depth_lm)  # a mask to get the components of the upper mantle
        mask_lm = (self.depths >= depth_lm)  # a mask to get the components of the lower mantle
        for i in range(N):
            Ediff = Ediffs[i]
            Edisl = Edisls[i]
            d = ds[i]
            if Ediff < 0.0 or Edisl < 0.0 or d <= 0.0:
                continue
            diffusion_creep['E'] = Ediff
            diffusion_creep['d'] = d
            dislocation_creep['E'] = Edisl
            dislocation_creep['d'] = d
            diffusion_creep['Edev'] = Ediff - diff_orig['E']
            dislocation_creep['Edev'] = Edisl - disl_orig['E']
            if (abs(diffusion_creep['Edev'] > Ediff_sigma) or abs(dislocation_creep['Edev'] > Edisl_sigma)):
                # check Edev is in range
                continue

            # 250km
            depth1 = 250e3
            eta250 = 5e19
            T1 = T_func(depth1)
            P1 = P_func(depth1)
            # diffusion creep
            # compute V, viscosities from both creeps are equal
            Vdiff = CreepComputeV(diffusion_creep, strain_rate, P1, T1, 2*eta250, d=d)
            diffusion_creep['V'] = Vdiff
            Vdisl = CreepComputeV(dislocation_creep, strain_rate, P1, T1, 2*eta250, d=d, use_effective_strain_rate=True)
            dislocation_creep['V'] = Vdisl
            if Vdiff < 0.0 or Vdisl < 0.0:
                continue
            diffusion_creep['Vdev'] = Vdiff - diff_orig['V']
            dislocation_creep['Vdev'] = Vdisl - disl_orig['V']

            # 660 km 
            # diffusion creep
            T660 = T_func(depth_lm)
            P660 = P_func(depth_lm)
            eta_diff660 = CreepRheology(diffusion_creep, strain_rate, P660, T660)
            # dislocation creep
            eta_disl660 = CreepRheology(dislocation_creep, strain_rate, P660, T660, use_effective_strain_rate=True)
            eta660 = ComputeComposite(eta_diff660, eta_disl660)
            
            # other constraints
            depth2 = 300e3
            Ttemp = T_func(depth2)
            Ptemp = P_func(depth2)
            # diffusion creep
            eta_diff300 = CreepRheology(diffusion_creep, strain_rate, Ptemp, Ttemp)
            # dislocation creep
            eta_disl300 = CreepRheology(dislocation_creep, strain_rate, Ptemp, Ttemp, use_effective_strain_rate=True)

            if include_lower_mantle is not None:
                # lower mantle rheology
                diff_lm = diffusion_creep.copy()
                diff_lm['V'] = LowerMantleV(diffusion_creep['E'], T_lm_mean, P_lm_mean, lm_grad_T, lm_grad_P)
                diff_lm['A'] = CreepComputeA(diff_lm, strain_rate, P660, T660, eta660*include_lower_mantle)

            
            # 1000km integral
            eta_diff = CreepRheology(diffusion_creep, strain_rate, self.pressures, self.temperatures)
            eta_disl = CreepRheology(dislocation_creep, strain_rate, self.pressures, self.temperatures, use_effective_strain_rate=True)
            eta = ComputeComposite(eta_diff, eta_disl)
            lith_depth = 100e3
            integral_depth = 1400e3
            mask_integral = (self.depths > lith_depth) & (self.depths < integral_depth)
            integral_cores = 4 * np.pi * (radius - self.depths)**2.0
            # upper mantle
            mask_um_integral = (mask_um & mask_integral)
            integral_um = np.trapz(eta[mask_um_integral] * integral_cores[mask_um_integral], self.depths[mask_um_integral])
            # lower mantle
            integral_lm = 0.0
            if include_lower_mantle is not None:
                mask_lm_integral = (mask_lm & mask_integral)
                eta_diff_lm = CreepRheology(diff_lm, strain_rate, self.pressures, self.temperatures)
                integral_lm = np.trapz(eta_diff_lm[mask_lm_integral] * integral_cores[mask_lm_integral], self.depths[mask_lm_integral])
            else:
                integral_lm = 4.0 / 3 * np.pi * ((radius - depth_lm)**3.0 - (radius - integral_depth)**3.0) * eta660 * 30.0 # assume 30 times jump
            volume = 4.0 / 3 * np.pi * ((radius - lith_depth)**3.0 - (radius - integral_depth)**3.0)
            average_eta = (integral_um + integral_lm) / volume

            # conditions:
            #   0: 300km
            #   1: 660km
            #   2: integral
            #   3: range of V
            eta660range = [4e20, 1e21]
            average_range = [0.65e21, 1.1e21]
            conds = [(eta_disl300 < eta_diff300),\
            (eta660 >= eta660range[0]) and (eta660 <= eta660range[1]),\
            (average_eta >= average_range[0]) and (average_eta <= average_range[1]),\
            abs(diffusion_creep['Vdev']) < Vdiff_sigma and abs(dislocation_creep['Vdev']) < Vdisl_sigma]
            # failed info
            # print(conds)
            condition_indexes = [0, 1, 3]
            cond_combined = True
            for i in condition_indexes:
                cond_combined = (cond_combined and conds[i])
            
            if cond_combined:
                constrained_ds.append(d)
                constrained_Vdisls.append(Vdisl)
                constrained_Vdiffs.append(Vdiff)
            else:
                constrained_rheologies.append({'diff': diffusion_creep.copy(), 'disl': dislocation_creep.copy(),\
                'average_upper_region': average_eta})


        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(constrained_Ediffs, constrained_Edisls, constrained_ds)
        ax.set_xlabel('Ediff [J/mol]')
        ax.set_ylabel('Edisl [J/mol]')
        ax.set_zlabel('d [um]')
        fig_name = 'constrained_rheology_%s_N%d.png' % (rheology, N)
        fig_path = os.path.join(fig_dir, fig_name)
        fig.savefig(fig_path)
        print("New figure: %s" % fig_path)
        plt.close()

        # plot profiles
        save_profile = kwargs.get('save_profile', 0)
        i = 0  # index
        for constrained_rheology in constrained_rheologies:
            # dump json
            Vdiff = constrained_rheology['diff']['V']
            Vdisl = constrained_rheology['disl']['V']
            Ediff = constrained_rheology['diff']['E']
            Edisl = constrained_rheology['disl']['E']
            d = constrained_rheology['disl']['d']
            json_path = os.path.join(fig_dir, "constrained_profile_Ediff%.4e_Edisl%.4e_Vdiff%.4e_Vdisl%.4e_d%.4e.json" % (Ediff, Edisl, Vdiff, Vdisl, d))
            with open(json_path, 'w') as fout:
                json.dump(constrained_rheology, fout)
                print("[%d / %d], New json: %s" % (i, len(constrained_rheologies), json_path))
            #  save profile
            if save_profile == 1:
                eta_diff = CreepRheology(constrained_rheology['diff'], strain_rate, self.pressures, self.temperatures)
                eta_disl = CreepRheology(constrained_rheology['disl'], strain_rate, self.pressures, self.temperatures, use_effective_strain_rate=True)
                eta = ComputeComposite(eta_diff, eta_disl)
                # plots
                fig, axs = plt.subplots(1, 2, figsize=(10, 5))
                color = 'tab:blue'
                axs[0].plot(self.pressures/1e9, self.depths/1e3, color=color, label='pressure')
                axs[0].set_ylabel('Depth [km]') 
                axs[0].set_xlabel('Pressure [GPa]', color=color) 
                # axs[0].invert_yaxis()
                if include_lower_mantle is None:
                    ylim=[660.0, 0.0]
                else:
                    ylim=[2890, 0.0]
                axs[0].set_ylim(ylim)
                # ax2: temperature
                color = 'tab:red'
                ax2 = axs[0].twiny()
                ax2.set_ylim(ylim)
                ax2.plot(self.temperatures, self.depths/1e3, color=color, label='temperature')
                ax2.set_xlabel('Temperature [K]', color=color) 
                # second: viscosity
                #   upper mantle
                axs[1].semilogx(eta_diff[mask_um], self.depths[mask_um]/1e3, 'c', label='diffusion creep')
                axs[1].semilogx(eta_disl[mask_um], self.depths[mask_um]/1e3, 'g', label='dislocation creep(%.2e)' % strain_rate)
                axs[1].semilogx(eta[mask_um], self.depths[mask_um]/1e3, 'r--', label='Composite')
                axs[1].set_xlim([1e18,1e25])
                axs[1].set_ylim(ylim)
                # axs[1].invert_yaxis()
                axs[1].grid()
                axs[1].set_ylabel('Depth [km]')
                axs[1].legend()
                if include_lower_mantle:
                    # include lower mantle info as title:
                    diff_lm = constrained_rheology['diff_lm']
                    eta_diff_lm = CreepRheology(constrained_rheology['diff_lm'], strain_rate, self.pressures, self.temperatures)
                    axs[1].semilogx(eta_diff_lm[mask_lm], self.depths[mask_lm]/1e3, 'c')
                    title_str = "lm_jump%.2d_Vlm%.4e_avg%.4e" % (include_lower_mantle, diff_lm['V'], constrained_rheology['average_upper_region'])
                    axs[1].set_title(title_str)
                # save figure
                fig_path = os.path.join(fig_dir, "constrained_profile_Ediff%.4e_Edisl%.4e_Vdiff%.4e_Vdisl%.4e_d%.4e.png" % (Ediff, Edisl, Vdiff, Vdisl, d))
                fig.savefig(fig_path)
                print("[%d / %d], New figure: %s" % (i, len(constrained_rheologies), fig_path))
                plt.close()
            i = i + 1




def GetRheology(rheology):
    '''
    read rheology parameters, and account for effects of water if it is a wet rheology
    '''
    RheologyPrm = RHEOLOGY_PRM()
    if not hasattr(RheologyPrm, rheology + "_diff") and hasattr(RheologyPrm, rheology + "_disl"):
        raise ValueError("RHEOLOGY_PRM object doesn't have attribute %s_diff" % rheology)
    diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
    dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
    try:
        _ = diffusion_creep['wet']
    except KeyError:
        pass
    else:
        ### effects of water accounted, see Magali's file explain_update_modHK03_rheology eq(5)
        water_creep = getattr(RheologyPrm, "water")
        diffusion_creep['A'] = diffusion_creep['A'] / (water_creep['A'] ** diffusion_creep['r'])
        diffusion_creep['V'] = diffusion_creep['V'] - water_creep['V'] * diffusion_creep['r']
        diffusion_creep['E'] = diffusion_creep['E'] - water_creep['E'] * diffusion_creep['r']
        dislocation_creep['A'] = dislocation_creep['A'] / (water_creep['A'] ** dislocation_creep['r'])
        dislocation_creep['V'] = dislocation_creep['V'] - water_creep['V'] * dislocation_creep['r']
        dislocation_creep['E'] = dislocation_creep['E'] - water_creep['E'] * dislocation_creep['r']
    return diffusion_creep, dislocation_creep


def Config(_kwargs, _name, _default):
    """    def ReadProfile(self, file_path):
        self.depths, self.pressures, self.temperatures = ReadAspectProfile(file_path)ble value and assign default if not found
    """
    try:
        value = _kwargs[_name]
    except KeyError:
        value = _default
    return value


def CreepStress(creep, strain_rate, P, T, d, Coh):
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
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # calculate B
    B = A * d**(-p) * Coh**r
    return (strain_rate / B)**(1.0 / n) * np.exp((E + P * V) / (n * R * T))


def CreepRheology(creep, strain_rate, P, T, d=None, Coh=None, **kwargs):
    """
    def CreepRheology(creep, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep['d']
    if Coh is None:
        Coh = creep['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    eta = 1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6

    return eta


def CreepRheology_v1(creep, strain_rate, P, T, d=None, Coh=None, **kwargs):
    """
    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Previously, there is a typo in the F factor
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+1)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep['d']
    if Coh is None:
        Coh = creep['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    eta = 1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6

    return eta


def CreepComputeV(creep, strain_rate, P, T, eta, d=None, Coh=None, **kwargs):
    """
    Calculate V based on other parameters 
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep['d']
    if Coh is None:
        Coh = creep['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    exponential = eta / (1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp(E / (n * R * T)) * 1e6)
    V = n * R * T * np.log(exponential) / P
    return V


def CreepComputeA(creep, strain_rate, P, T, eta, d=None, Coh=None, **kwargs):
    """
    def CreepRheology(creep, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep['d']
    if Coh is None:
        Coh = creep['Coh']
    # calculate B
    B = (F/eta)**n * strain_rate**(1-n) * np.exp((E+P*V)/(R*T)) * (1e6)**n
    A = B * d**p * Coh**(-r)
    return A


def CreepComputeA_v1(creep, strain_rate, P, T, eta, d=None, Coh=None, **kwargs):
    """
    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    Here I tried to input the right value for the F factor
    """
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+1)/2/n))
    else:
        F = 1.0

    if d is None:
        d = creep['d']
    if Coh is None:
        Coh = creep['Coh']
    # calculate B
    B = (F/eta)**n * strain_rate**(1-n) * np.exp((E+P*V)/(R*T)) * (1e6)**n
    A = B * d**p * Coh**(-r)
    return A


def CreepRheologyInAspectViscoPlastic(creep, strain_rate, P, T):
    """
    def CreepRheologyInAspectVisoPlastic(creep, strain_rate, P, T)

    Calculate viscosity by way of Visco Plastic module in aspect
    flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: m
     - Return value: Pa*s
    """
    A = creep['A']
    m = creep['m']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    d = creep['d']
    # calculate B
    return 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))


def Convert2AspectInput(creep, **kwargs):
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
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    d = creep['d']
    Coh = creep['Coh']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+2)/2/n))
    else:
        F = 1.0
    # prepare values for aspect
    aspect_creep = {}
    # stress in the original equation is in Mpa, grain size is in um
    aspect_creep['A'] = 1e6**(-p) * (1e6)**(-n) * Coh**r * A / F**n  # F term: use effective strain rate
    aspect_creep['d'] = d / 1e6
    aspect_creep['n'] = n
    aspect_creep['m'] = p
    aspect_creep['E'] = E
    aspect_creep['V'] = V
    return aspect_creep


def Convert2AspectInput_v1(creep, **kwargs):
    """
    Viscosity is calculated by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6
    while in aspect, flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    In this version, I am trying to take care of the F factor correctly
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
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    d = creep['d']
    Coh = creep['Coh']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+1)/2/n))
    else:
        F = 1.0
    # prepare values for aspect
    aspect_creep = {}
    # stress in the original equation is in Mpa, grain size is in um
    aspect_creep['A'] = 1e6**(-p) * (1e6)**(-n) * Coh**r * A / F**n  # F term: use effective strain rate
    aspect_creep['d'] = d / 1e6
    aspect_creep['n'] = n
    aspect_creep['m'] = p
    aspect_creep['E'] = E
    aspect_creep['V'] = V
    return aspect_creep


def ConvertFromAspectInput(aspect_creep, **kwargs):
    """
    Viscosity is calculated by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6
    while in aspect, flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T)).
    Here I convert backward from the flow law used in aspect
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
    A = aspect_creep['A']
    m = aspect_creep['m']
    n = aspect_creep['n']
    E = aspect_creep['E']
    V = aspect_creep['V']
    d = aspect_creep['d']
    # compute value of F(pre factor)
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 1 / (2**((n-1)/n)*3**((n+1)/2/n))
    else:
        F = 1.0
    # prepare values for aspect
    creep = {}
    # stress in the original equation is in Mpa, grain size is in um
    creep['A'] = 1e6**m * 1e6**n * A * F**n  # F term: use effective strain rate
    creep['d'] = d * 1e6
    creep['n'] = n
    creep['p'] = m
    creep['E'] = E
    creep['V'] = V
    creep['r'] = 0.0  # assume this is not dependent on Coh
    creep['Coh'] = 1000.0
    return creep


def GetLowerMantleRheology(upper_mantle_creep_method, jump, T, P, **kwargs):
    """
    get flow law parameters in lower mantle based on upper mantle viscosity and jump in viscosity
    variables:
     - jump: viscosity jump at 660km
     - T: temperature at 660km
     - P: pressure at 660km
    """
    # extra inputs
    strategy = kwargs.get('strategy', 'A')
    V1 = kwargs.get('V1', upper_mantle_creep_method['V'])
    # read upper mantle values
    A = upper_mantle_creep_method['A']
    m = upper_mantle_creep_method['m']
    n = upper_mantle_creep_method['n']
    E = upper_mantle_creep_method['E']
    V = upper_mantle_creep_method['V']
    d = upper_mantle_creep_method['d']

    lower_mantle_creep_method = dict(upper_mantle_creep_method)
    lower_mantle_creep_method['V'] = V1
    if strategy == 'A':
        lower_mantle_creep_method['A'] = jump**(-n) * A * math.exp(P * (V1 - V) / (R * T))
    elif strategy == 'composite':
        # composite: prescribe a P and a V here, as well as pressure and temperature at 660 km depth. 
        # The composite viscosity of upper mantle is used as the base,
        # and a upper_lower_viscosity factor will be multiplied on that.
        eta_lower = kwargs['eta'] * jump # viscosity at 660 km
        lower_mantle_creep_method['V'] = V1
        lower_mantle_creep_method['A'] = 0.5/eta_lower * d**m * math.exp((E + P*V1) / (R*T)) # diffusion, thus n = 1
    elif strategy == 'c12':
        # c12: use P and V value from cizcova et al 2012, compute A using value of pressure and temperature at 660 km depth
        eta_lower = kwargs['eta'] * jump # viscosity at 660 km
        V1 = 1.1e-6
        E = 2e5
        lower_mantle_creep_method['V'] = V1
        lower_mantle_creep_method['E'] = E
        lower_mantle_creep_method['A'] = 0.5/eta_lower * d**m * math.exp((E + P*V1) / (R*T)) # diffusion, thus n = 1
    elif strategy == 'c12_const':
        # c12_const: only use the value of constrainted, 3-4e22
        eta_lower = 3.5e22
        lower_mantle_creep_method['V'] = 0.0
        lower_mantle_creep_method['E'] = 0.0
        lower_mantle_creep_method['m'] = 0.0
        lower_mantle_creep_method['A'] = 0.5/eta_lower
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
    time_step = 0
    i0 = DepthAverage.time_step_indexes[time_step][-1] * DepthAverage.time_step_length
    if time_step == len(DepthAverage.time_step_times) - 1:
        # this is the last step
        i1 = DepthAverage.data.shape[0]
    else:
        i1 = DepthAverage.time_step_indexes[time_step + 1][0] * DepthAverage.time_step_length
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
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    strain_rate = kwargs.get('strain_rate', 1e-15)
    # grain size
    try:
        d = kwargs['d']
        diffusion_creep['d'] = d
        dislocation_creep['d'] = d
    except KeyError:
        d = diffusion_creep['d']

    
    # diffusion creep
    implementation = kwargs.get('implementation', 'LHY')
    if implementation == 'LHY':
        eta_diff = CreepRheology(diffusion_creep, strain_rate, pressures, temperatures)
        eta_annotation += rheology
    elif implementation == 'MB':
        # use magali's implementation
        coh = 1000
        water = 'wet' # 'wet, dry, con'
        mod = 'new'  # orig, new
        Edev = 'mid'
        Vdev = 'mid'
        eta_diff = visc_diff_HK(temperatures,pressures,d,coh,water,mod,Edev,Vdev)
        eta_annotation += '%s_%s_E%s_V%s' % (water,mod, Edev, Vdev)
    else:
        raise CheckValueError('%s is not a valid implementation' % implementation)
    # dislocation creep
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
    

def LowerMantleV(E, Tmean, Pmean, grad_T, grad_P):
    '''
    compute the value of activation volume for the lower mantle
    based on the criteria of a nearly constant viscosity
    '''    
    V = E * grad_T / (grad_P * Tmean - Pmean * grad_T)
    return V


def ConstrainASPECT(file_path, **kwargs):
    '''
    Figure out contrain of rheology by random walk
    Inputs:
        file_path(str): a profile from ASPECT
    '''
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(file_path)
    # do a random walk
    save_profile = kwargs.get('save_profile', 0)
    include_lower_mantle = kwargs.get('include_lower_mantle', None)
    version = kwargs.get('version', 0)
    if version == 0:
        Operator.ConstrainRheology_v0(save_profile=save_profile, include_lower_mantle=include_lower_mantle)
    elif version == 1:
        Operator.ConstrainRheology_v1(save_profile=save_profile, include_lower_mantle=include_lower_mantle)


def DeriveMantleRheology(file_path, **kwargs):
    '''
    Derive a Mantle rheology profile following certain procedures
    Inputs:
        file_path(str): a profile from ASPECT
    '''
    mantle_rheology_scheme = kwargs.get("rheology", "HK03_wet_mod")
    diff = kwargs.get("diff", [1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    assert(type(diff) == list)
    dAdiff_ratio = diff[0]
    dEdiff = diff[1]
    dVdiff = diff[2]
    dAdisl_ratio = diff[3]
    dEdisl = diff[4]
    dVdisl = diff[5]
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(file_path)
    # do a random walk
    save_profile = kwargs.get('save_profile', 0)
    include_lower_mantle = kwargs.get('include_lower_mantle', None)
    version = kwargs.get('version', 0)
    if version == 0:
        Operator.MantleRheology_v0(rheology=mantle_rheology_scheme,save_profile=save_profile)
    elif version == 1:
        Operator.MantleRheology_v1(rheology=mantle_rheology_scheme,save_profile=save_profile,\
            dEdiff=dEdiff, dEdisl=dEdisl, dVdiff=dVdiff, dVdisl=dVdisl,\
                dAdiff_ratio=dAdiff_ratio, dAdisl_ratio=dAdisl_ratio)
    else:
        raise CheckValueError('%d is not a valid version of Mantle Rheology' % version)

###
# functions for deriving the strenght profile
###
class STRENGTH_PROFILE(RHEOLOGY_OPR):

    def __init__(self):
        RHEOLOGY_OPR.__init__(self)
        self.Sigs = None
        self.Sigs_viscous = None
        self.Zs = None
        self.Etas = None
        self.Computed = False
        self.creep_type = None
        self.creep = None

    def Execute(self, **kwargs):
        '''
        Compute the strength profile
        '''
        year = 365 * 24 * 3600.0
        self.creep_type = kwargs.get('creep_type', 'diff')
        self.creep = getattr(self, self.creep_type)
        plastic = self.plastic
        assert(self.creep != None)
        assert(plastic != None)
        averaging = kwargs.get('averaging', 'harmonic')
        strain_rate = kwargs.get('strain_rate', 1e-14)
        # rheology_prm = RHEOLOGY_PRM()
        # self.plastic = rheology_prm.ARCAY17_plastic
        # self.dislocation_creep = rheology_prm.ARCAY17_disl
        Zs = np.linspace(0.0, 40e3, 100)
        Tliths = temperature_halfspace(Zs, 40e6*year, Tm=1573.0) # adiabatic temperature
        Ts = 713 * Zs / 78.245e3  + 273.14# geotherm from Arcay 2017 pepi, figure 3d 2
        Ps = pressure_from_lithostatic(Zs, Tliths)
        # plastic self.plastic
        if plastic["type"] == "stress dependent":
            Sigs_plastic = StressDependentYielding(Ps, plastic["cohesion"], plastic["friction"], plastic["ref strain rate"], plastic["n"], strain_rate)
        elif plastic["type"] == "Coulumb":
            Sigs_plastic = CoulumbYielding(Ps, plastic["cohesion"], plastic["friction"])
        eta_plastic = Sigs_plastic / 2.0 / strain_rate
        # viscous stress
        Sigs_viscous = CreepStress(self.creep, strain_rate, Ps, Ts, 1e4, 1000)
        eta_viscous = Sigs_viscous / 2.0 / strain_rate
        Sigs = np.minimum(Sigs_plastic, Sigs_viscous)
        if averaging == 'harmonic':
            Etas = eta_plastic * eta_viscous / (eta_viscous + eta_plastic)
        else:
            raise NotImplementedError()
        self.Sigs = Sigs
        self.Sigs_viscous = Sigs_viscous
        self.Etas = Etas
        self.Zs = Zs
        self.computed = True

    def PlotStress(self, **kwargs):
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        label_viscous = kwargs.get('label_viscous', None)
        _color = kwargs.get('color', 'b')
        if ax == None:
            raise NotImplementedError()
        # make plots
        ax.plot(self.Sigs/1e6, self.Zs/1e3, color=_color, label=label)
        ax.plot(self.Sigs_viscous/1e6, self.Zs/1e3, '--', color=_color, label=label_viscous)
        ax.set_xlim([0.0, 100.0])
        ax.set_xlabel("Second invariant of the stress tensor (MPa)")
        ax.set_ylabel("Depth (km)")
    
    def PlotViscosity(self, **kwargs):
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        _color = kwargs.get('color', 'b')
        if ax == None:
            raise NotImplementedError()
        # plot viscosity
        ax.semilogx(self.Etas, self.Zs/1e3, color=_color, label=label)
        ax.set_xlabel("Viscosity (Pa * s)")
        ax.set_ylabel("Depth (km)")
        

def PlotShearZoneStrengh(Operator, fig_path_base):
    '''
    Plot the shear zone strenght profile
    '''
    fig = plt.figure(tight_layout=True, figsize=[5, 10])
    gs = gridspec.GridSpec(2, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    strain_rates = [1e-13, 1e-14, 1e-15]
    colors = ['b', 'g', 'r']
    # 1e-13 
    i = 0
    for strain_rate in strain_rates:
        _color = colors[i]
        Operator.Execute(creep_type='disl', strain_rate=strain_rate)
        # plot stress
        label = "Strain Rate = %.1e" % strain_rate
        Operator.PlotStress(ax=ax0, color=_color, label_viscous=label)
        # plot viscosity
        Operator.PlotViscosity(ax=ax1, color=_color)
        i += 1
    ax0.invert_yaxis()
    ax0.legend()
    ax1.invert_yaxis()
    ax1.legend()
    # figure path
    fig_path = fig_path_base.split('.')[0]
    if Operator.creep_type == 'diff':
        fig_path += ("_diff_" + Operator.diff_type)
    elif Operator.creep_type == 'disl':
        fig_path += ("_disl_" + Operator.disl_type)
    elif Operator.creep_type == 'comp':
        fig_path += (("_diff_" + Operator.diff_type) + ("_disl_" + Operator.disl_type))
    else:
        raise TypeError("Wrong type of creep type %s" % Operator.creep_type)
    fig_path += '.' + fig_path_base.split('.')[1]
    fig.savefig(fig_path)
    print("figure saved: ", fig_path)


# yielding criteria
def Byerlee(P):
    '''
    byerlee's law for yielding
    Inputs:
        P (pressure, Pa) - lithostatic pressure
    '''
    if type(P) == float:
        sigma_n = P  # future: pore pressure
        if P < 200e6:
            tau = 0.85 * P
        else:
            tau = 0.6 * P + 60e6
    elif type(P) == np.ndarray:
        tau = np.zeros(P.shape)
        mask = (P < 200e6)
        tau[mask] = 0.85 * P[mask]
        tau[~mask] = 0.6 * P[~mask] + 60e6
    else:
        raise TypeError("Wrong type of entry")
    return tau


def StressDependentYielding(P, cohesion, friction, strain_rate_ref, n, strain_rate):
    '''
    a yielding criteria that include a stress dependence on strain rate
    '''
    tau_y = cohesion + friction * P
    tau = tau_y * (strain_rate / strain_rate_ref) ** (1.0/n)
    return tau

def CoulumbYielding(P, cohesion, friction):
    '''
    a yielding criteria that include a stress dependence on strain rate
    '''
    tau = cohesion + friction * P
    return tau


def pressure_from_lithostatic(z,Tad):
    '''
    A lithostatic pressure profile
    Inputs:
	z (float) - depth in m
	Tad (float) - adiabatic temperature
    '''
    # Density Profile
    refrho = 3300  # kg/m^3
    refT = 1673        # K
    alpha = 3.1e-5  # 1/K
    g = 9.81 # m/s^2
    density = refrho*(1-alpha*(Tad-refT))
    # start loop at 1 because P[0] = 0
    dz = z[1]-z[0]
    P = np.zeros(np.size(z))
    for i in range(1, np.size(z)):
        P[i] = P[i-1] + 0.5*(density[i]+density[i-1])*g*dz
    return P

def temperature_halfspace(z, t, **kwargs):
    '''
    temperature from a half-space cooling model
    Inputs:
	z (float) - depth (m)
	t - age (s)
    kwargs (dict):
        Tm (float) - mantle temperature
    '''
    # Physical constants
    kappa = 1e-6  # thermal diffusivity (m^2/s)
    T_s = 273  # surface temperature (K)
    T_m = kwargs.get("Tm", 1673) # mantle temperature (K)
    T = T_s + (T_m - T_s)*erf(z/(2*np.sqrt(kappa*t)))
    return T




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
                        default=None,
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
    parser.add_argument('-sf', '--save_profile', type=int,
                        default=1,
                        help='Save profile when constraining rheology')
    parser.add_argument('-ilm', '--include_lower_mantle', type=float,
                        default=None,
                        help='Include lower mantle in the computation of rheology')
    parser.add_argument('-v', '--version', type=int,
                        default=0,
                        help='Version')
    parser.add_argument('-df', '--diff', type=str,
                        default="1.0,0.0,0.0,1.0,0.0,0.0",
                        help='Differences from the default values in a rheology')
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
        diffusion_creep, dislocation_creep = GetRheology(rheology)
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
        fig_path = os.path.join(RESULT_DIR, 'along_profile_rhoelogy.png')
        depths, pressures, temperatures = ReadAspectProfile(arg.inputs)
        if arg.grain_size is not None:
            PlotAlongProfile(depths, pressures, temperatures, fig_path, rheology=arg.rheology, d=arg.grain_size, implementation=arg.implementation)
        else:
            # version without a grain size, use grain size predifined in the rheology
            PlotAlongProfile(depths, pressures, temperatures, fig_path, rheology=arg.rheology, implementation=arg.implementation)
    
    elif _commend == 'plot_along_aspect_profile_with_json':
        fig_path = os.path.join(RESULT_DIR, 'along_profile_rhoelogy.png')
        depths, pressures, temperatures = ReadAspectProfile(arg.inputs)
        PlotAlongProfileJson(depths, pressures, temperatures, arg.json, fig_path)

    elif _commend == 'constrain_aspect_rheology':
        ConstrainASPECT(arg.inputs, save_profile=arg.save_profile, include_lower_mantle=arg.include_lower_mantle, version=arg.version)

    elif _commend == 'derive_mantle_rheology':
        diff_str = arg.diff.split(',')
        diff = [float(entry) for entry in diff_str]
        DeriveMantleRheology(arg.inputs, save_profile=arg.save_profile, version=arg.version, rheology=arg.rheology, diff=diff)

    elif _commend == 'plot_strength_profile':
        Operator = STRENGTH_PROFILE()
        Operator.SetRheologyByName(disl='ARCAY17', plastic='ARCAY17')
        Operator.Execute(creep_type='disl')
        Sigs = Operator.Sigs
        Zs = Operator.Zs
        fig_path = os.path.join(ASPECT_LAB_DIR, "results", "strength_profile.png")
        PlotShearZoneStrengh(Operator, fig_path)
    else:
        raise CheckValueError('%s is not a valid commend' % _commend)



# run script
if __name__ == '__main__':
    main()


