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

  - Plot the shear-zone strength

        - plot for a specific strain rate
        python -m shilofue.Rheology plot_strength_profile -S 1e-15

        - plot it for a range of strain rate [1e-13, 1e-14, 1e-15] 
        python -m shilofue.Rheology plot_strength_profile -S -1.0

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
from matplotlib import gridspec, cm
from matplotlib import patches as mpatches 
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
from shilofue.FlowLaws import visc_diff_HK
from shilofue.ParsePrm import ParseFromDealiiInput, UpperMantleRheologyViscoPlastic, ReplacePhaseOption
from shilofue.ThermalModel import MANTLE_ADIABAT
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from shutil import rmtree, copy

R = 8.314

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')

# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

class CheckValueError(Exception):
    pass


class RHEOLOGY_PRM():
    """
    class for rheologies
    components and units:
        A (the prefactor) - MPa^(-n-r)*um**p/s
        n (stress dependence) - 1
        p (grain size dependence) - 1
        r (power of C_{OH}) - 1
        E (activation energy) - J / mol
        V (activation volume) - m^3 / mol
    Notes on the choice of the units:
        The unit of E and V are much easier to convert to UI.
        But for A, this code will handle the convertion, so the
        user only need to take the value for published flow laws.
    """
    def __init__(self):
        '''
        Initiation, initiate rheology parameters
        '''
        self.HK03_dry_disl = \
            {
                "A": 1.1e5,
                "p": 0.0,
                "r": 0.0,
                "n": 3.5,
                "E": 530e3,
                "V": 12e-6,
            }
        
        # dry diffusion creep in Hirth & Kohlstedt 2003)
        # note the V (activation energy) value has a large variation, here I
        # picked up a value the same as the wet value.
        self.HK03_dry_diff = \
            {
                "A": 1.5e9,
                "p": 3.0,
                "r": 0.0,
                "n": 1.0,
                "E": 375e3,
                "V": 4e-6,
            }

        # dislocation creep in Hirth & Kohlstedt 2003
        # with constant Coh
        self.HK03_disl = \
            {
                "A": 90,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
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
            }
        
        # dislocation creep in Hirth & Kohlstedt 2003
        # with varied Coh
        self.HK03_w_disl = \
            {
                "A": 1600.0,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 520e3,
                "V": 22e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03_w_diff = \
            {
                "A" : 2.5e7,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3,
                "V" : 10e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0
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
                "A" : 2.85e5,  # note: their number in the 2017 appendix is wrong,
                "p" : 3.0, #  but it's right in the 2016 paper.
                "r" : 1.0,
                "n" : 1.0,
                "E" : 317e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        # 'wet' indicates this has to applied with a rheology of water
        self.HK03_wet_mod_diff = \
            {
                "A" : 7.1768184e6,  
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
                "A" : 457.142857143,
                "p" : 0.0,
                "r" : 1.2,
                "n" : 3.5,
                "E" : 520e3,
                "V" : 24e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet" : 1.0
            }
        
        # modified creep laws from Hirth & Kohlstedt 2003
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        # 'wet' indicates this has to applied with a rheology of water
        self.HK03_wet_mod_ln_diff = \
            {
                "A" : 7.1768184e6,  
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 375e3,
                "V" : 23e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod_ln_disl = \
            {
                "A" : 457.142857143,
                "p" : 0.0,
                "r" : 1.2,
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
        # for detail, refer to magali's explain_update_modHK03_rheology.pdf file
        # 'wet' indicates this has to applied with a rheology of water
        # In the version, I modified the value of r, compared to the first version
        self.HK03_wet_mod2_diff = \
            {
                # "A" : 10**6.9,  # MPa^(-n-r)*um**p/s
                "A" : 7.1768e6,  # MPa^(-n-r)*um**p/s
                "p" : 3.0,
                "r" : 0.8, # 1.0 -> 0.8
                "n" : 1.0,
                "E" : 375e3,
                "V" : 23e-6,
                "d" : 1e4,
                "Coh" : 1000.0,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }

        self.HK03_wet_mod2_disl = \
            {
                "A" : 10**2.65,
                "p" : 0.0,
                "r" : 1.2,  # 1.0 -> 1.2
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
        # this is specifically the one I used for the TwoD models.
        # Combined with the usage of function "MantleRheology", then same rheology
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
        self.ARCAY17_brittle = \
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

        # todo_mineral
        # Basalt rheology from Shelton and Tullis 1981 
        # and Hacker and Christie 1990
        # the diffusion creep of this rheology is missing
        # in literatures
        self.ST1981_basalt_diff = None

        self.ST1981_basalt_disl = \
            {
                "A" : 1.0e-4,
                "p" : 0.0,
                "r" : 0.0,  # not dependent on the "Coh"
                "n" : 3.5,
                "E" : 250e3,
                "V" : 0.0,
                "d" : 1e4, # not dependent on d
            }

        # Quartz rheology from Ranalli and Murphy 1987.
        # the diffusion creep of this rheology is missing
        # in literatures
        self.ST1981_basalt_diff = None

        self.RM1987_quartz_disl = \
            {
                "A" : 6.8e-6,
                "p" : 0.0,
                "r" : 0.0,  # not dependent on the "Coh"
                "n" : 3,
                "E" : 156e3,
                "V" : 0.0
            }
        
        self.KK1987_quartz_disl = \
            {
                "A" : 3.2e-4,
                "p" : 0.0,
                "r" : 0.0,  # not dependent on the "Coh"
                "n" : 2.3,
                "E" : 154e3,
                "V" : 8e-6
            }
        
        self.Ranali_95_anorthite_75_diff = None
        
        self.Ranali_95_anorthite_75_disl = \
            {
                "A" : 3.3e-4,
                "p" : 0.0,
                "r" : 0.0,  # not dependent on the "Coh"
                "n" : 3.2,
                "E" : 238e3,
                "V" : 8e-6
            }

        self.Rybachi_06_anorthite_wet_diff = \
            {
                "A" : 0.2,  # note, 10^(-0.7), less than 1 digit accuracy, as 10^(0.1) = 1.25
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 159e3,
                "V" : 38e-6,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }
        
        self.Rybachi_06_anorthite_wet_disl = \
            {
                "A" : 1.6,  # note, 10^(0.2), less than 1 digit accuracy, as 10^(0.1) = 1.25
                "p" : 0.0,
                "r" : 1.0,
                "n" : 3.0,
                "E" : 345e3,
                "V" : 38e-6,
                "wet": 1.0  # I use this to mark this is a wet rheology, so I need to account for V and E for water later.
            }
        
        self.Rybachi_06_anorthite_dry_diff = \
            {
                "A" : 1.26e12,  # note, 10^(12.1), less than 1 digit accuracy, as 10^(0.1) = 1.25
                "p" : 3.0,
                "r" : 0.0,  # dry, not dependent on fugacity
                "n" : 1.0,
                "E" : 460e3,
                "V" : 24e-6
            }
        
        self.Rybachi_06_anorthite_dry_disl = \
            {
                "A" : 5.01e12,  # note, 10^(12.7), less than 1 digit accuracy, as 10^(0.1) = 1.25
                "p" : 0.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 3.0,
                "E" : 641e3,
                "V" : 24e-6
            }
        
        self.Dimanov_Dresen_An50Di35D_wet_diff = \
        {
                # diffusion creep for a 35 mu*m grain size
                "A" : 5488000000.0,  #   1.28e-1 * (1e6) / (35)^(-3)
                "p" : 3.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 1.0,
                "E" : 316e3,
                "V" : 0.0  # not present in the table
        }
        
        self.Dimanov_Dresen_An50Di35D_wet_disl = \
        {
            # this is actually the An50DiD in table 3b
            # since the dislocation creep is not grain size sensitive
                "A" : 10174679.0993,  #  / 1.54e-17 * (1e6)^3.97
                "p" : 0.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 3.97,
                "E" : 556e3,
                "V" : 0.0  # not present in the table
        }

        self.Dimanov_Dresen_An50Di35D_dry_diff = \
        {
                # diffusion creep for a 35 mu*m grain size
                "A" : 5.1879e13,  #  1.21e3 * (1e6) / (35)^(-3)
                "p" : 3.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 1.0,
                "E" : 436e3,
                "V" : 0.0  # not present in the table
        }
        
        self.Dimanov_Dresen_An50Di35D_dry_disl = \
        {
            # this is actually the An50DiD in table 3b
            # since the dislocation creep is not grain size sensitive
                "A" : 8.1840692e+12,  #  / 2.71e-12 * (1e6)^4.08
                "p" : 0.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 4.08,
                "E" : 723e3,
                "V" : 0.0  # not present in the table
        }
        
        
        self.Dimanov_Dresen_An50Di45D_dry_diff = \
        {
                # diffusion creep for a 45 mu*m grain size
                "A" : 6.187e15,  # 6.79e4 * 1e6 / (45)^(-3.0)
                "p" : 3.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 1.0,
                "E" : 496e3,
                "V" : 0.0  # not present in the table
        }
        
        self.Dimanov_Dresen_An50Di45D_dry_disl = \
        {
            # this is actually the An50DiD in table 3b
            # since the dislocation creep is not grain size sensitive
                "A" : 8.1840692e+12,  #  / 2.71e-12 * (1e6)^4.08
                "p" : 0.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 4.08,
                "E" : 723e3,
                "V" : 0.0  # not present in the table
        }

        self.Rybachi_2000_An100_dry_diff = \
        {
                # diffusion creep, for the An100 in table 3
                # note that this is marked as dry, but there
                # is 640 ppm H/Si in the synthetic anorthite
                "A" : 1.258925e12, # 10^12.1
                "p" : 3.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 1.0,
                "E" : 467e3,
                "V" : 0.0  # not present in the table
        }
        
        self.Rybachi_2000_An100_dry_disl = \
        {
            # dislocation creep, for the An100 in table 3
                "A" : 5.01187e12,  # 10^12.7
                "p" : 0.0,
                "r" : 0.0, # dry, not dependent on fugacity
                "n" : 3.0,
                "E" : 648e3,
                "V" : 0.0  # not present in the table
        }

        # todo_peierls 
        self.MK10_peierls = \
        {
            'q': 1.0,
            'p': 0.5,
            'n': 2.0,
            'sigp0': 5.9e3,    				# MPa (+/- 0.2e3 Pa)
            'A': 1.4e-7,      # s^-1 MPa^-2
            'E': 320e3,      				# J/mol (+/-50e3 J/mol)
            'V' : 0.0  # not dependent on the pressure
        }
        
        self.Idrissi16_peierls = \
        {
            'q': 2.0,
            'p': 0.5,
            'n': 0.0,
            'sigp0': 3.8e3,    				# MPa (+/- 0.2e3 Pa)
            'A': 1e6,      # s^-1 MPa^-2
            'E': 566e3,      				# J/mol (+/-50e3 J/mol)
            'V' : 0.0  # not dependent on the pressure
        }


        self.Byerlee_brittle = \
        {
            "type": "Byerlee"
        }

    def get_rheology(self, _name, _type):
        '''
        read rheology parameters, and account for effects of water if it is a wet rheology
        '''
        assert(_type in ['diff', 'disl', 'brittle'])
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

class RHEOLOGY_OPT(Utilities.JSON_OPT):
    '''
    Define a complex class for using the json files for testing the rheologies.
    '''
    def __init__(self):
        '''
        initiation
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Type of diffusion creep", str, ["diffusion"], "", nick='diffusion')
        self.add_key("Type of dislocation creep", str, ["dislocation"], "", nick='dislocation')
        self.add_key("Differences in ratio of the prefactor for diffusion creep", float, ["diffusion prefactor difference ratio"], 1.0, nick='dA_diff_ratio')
        self.add_key("Differences of the activation energy for diffusion creep", float, ["diffusion activation energy difference"], 0.0, nick='dE_diff')
        self.add_key("Differences of the activation volume for diffusion creep", float, ["diffusion activation volume difference"], 0.0, nick='dV_diff')
        self.add_key("Differences in ratio of the prefactor for dislocation creep", float, ["dislocation prefactor difference ratio"], 1.0, nick='dA_disl_ratio')
        self.add_key("Differences of the activation energy for dislocation creep", float, ["dislocation activation energy difference"], 0.0, nick='dE_disl')
        self.add_key("Differences of the activation volume for dislocation creep", float, ["dislocation activation volume difference"], 0.0, nick='dV_disl')
        # todo_r_json
        self.add_key("Grain size in mu m", float, ["grain size"], 10000.0, nick='d')
        self.add_key("Coh in /10^6 Si", float, ["coh"], 1000.0, nick='coh')
        self.add_key("fh2o in MPa", float, ["fh2o"], -1.0, nick='fh2o')
    
    def check(self):
        '''
        check values are validate
        '''
        RheologyPrm = RHEOLOGY_PRM()
        diffusion = self.values[0]
        if diffusion != "":
            diffusion_creep_key = (diffusion + "_diff") 
            assert(hasattr(RheologyPrm, diffusion_creep_key))
        dislocation = self.values[1]
        if dislocation != "":
            dislocation_creep_key = (dislocation + "_disl") 
            assert(hasattr(RheologyPrm, dislocation_creep_key))

    def to_RheologyInputs(self):
        '''
        '''

        diffusion = self.values[0]
        dislocation = self.values[1]
        dA_diff_ratio = self.values[2]
        dE_diff = self.values[3]
        dV_diff = self.values[4]
        dA_disl_ratio = self.values[5]
        dE_disl = self.values[6]
        dV_disl = self.values[7]
        d = self.values[8]
        coh = self.values[9]
        fh2o = self.values[10]
        use_coh = True
        if fh2o > 0.0:
            use_coh = False
        return diffusion, dislocation, dA_diff_ratio, dE_diff, dV_diff,\
        dA_disl_ratio, dE_disl, dV_disl, d, coh, fh2o, use_coh


class RHEOLOGY_PLOT_OPT(Utilities.JSON_OPT):
    '''

    '''
    def __init__(self):
        '''
        initiation
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Strain rate array", list, ["strain rates"], [1e-15], nick='strain_rates')
        self.add_key("Grain size array", list, ["grain sizes"], [1e4], nick='grain_sizes')
        self.add_key("type of plot", str, ["type"], "strain rate vs stress", nick='_type')

        # todo_r_json
        self.add_key("temperature (C)", float, ["T"], 1600.0, nick='temperature')


    def check(self):
        '''
        check
        '''
        _type = self.values[2]
        assert(_type in ["strain rate vs stress", "viscosity vs temperature"])
    
    def GetStrainRates(self):
        '''
        return the value of grain sizes
        '''
        return self.values[0]
    
    def GetGrainSizes(self):
        '''
        return the value of grain sizes
        '''
        return self.values[1]
    
    def GetType(self):
        '''
        return the type of plot
        '''
        return self.values[2]
    
    def GetT(self):
        '''
        return the temperature of plot
        '''
        return self.values[3]
    
        

class RHEOLOGY_JSON(Utilities.JSON_OPT):
    '''
    Define a complex class for using the json files for testing the rheologies.
    '''
    def __init__(self):
        '''
        initiation
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_features('Rheology options', ['rheologies'], RHEOLOGY_OPT)
        self.add_features('Plot options', ['plots'], RHEOLOGY_PLOT_OPT)
        self.add_key("path of plot", str, ["figure path"], "./foo.png", nick='figure_path')

    def check(self):
        rheology_features = self.values[0]
        plot_features = self.values[1]
        for plotOpt in plot_features:
            assert (len(plotOpt.GetGrainSizes()) in [1, len(rheology_features)])
            assert (len(plotOpt.GetStrainRates()) in [1, len(rheology_features)])
        
        figure_path = Utilities.var_subs(self.values[2])
        assert(os.path.isdir(os.path.dirname(figure_path)))

    def GetRheologyFeatures(self):
        '''
        Return the rheology features in an array
        '''
        rheology_features = self.values[0]
        return rheology_features
    
    def GetPlotFeatures(self):
        '''
        return the features of plot
        '''
        plot_features = self.values[1]
        return plot_features
    
    def GetFigurePath(self):
        '''
        return the path of the figure
        ''' 
        figure_path = Utilities.var_subs(self.values[2])
        return figure_path


class RHEOLOGY_OPR():
    '''
    rheology operation, do some complex staff
    Attributes:
        RheologyPrm: an initiation of the class RHEOLOGY_PRM
        depths(ndarray), pressures, temperatures: depth, pressure, temperature profile
            (all these 3 profiles are loaded from a depth_average outputs of ASPECT)
        peierls_type: a string, type of the peierls rheology
        peierls: the dictionary of variables to be used for the peierls rheology
    '''
    def __init__(self):
        '''
        Initiation
        '''
        self.RheologyPrm = RHEOLOGY_PRM()
        # set of variables for mantle profile
        self.depths = None
        self.pressures = None
        self.tempertures = None
        self.output_profile = None # for the figure plotted
        self.output_json = None
        self.output_aspect_json = None
        self.diff_type = None
        self.diff = None
        self.disl_type = None
        self.disl = None
        self.brittle_type = None
        self.brittle = None
        self.peierls_type = None
        self.peierls = None
        pass

    # todo_peierls
    def SetRheology(self, **kwargs):
        '''
        set rheology type with instances of rheology (i.e. a dictionary)
        '''
        self.diff = kwargs.get('diff', None)
        self.disl = kwargs.get('disl', None)
        self.brittle = kwargs.get('brittle', None)
        self.peierls = kwargs.get('peierls', None)
        pass
    
    def SetRheologyByName(self, **kwargs):
        '''
        set rheology type with instances of rheology (i.e. a dictionary)
        '''
        diff_type = kwargs.get('diff', None)
        disl_type = kwargs.get('disl', None)
        brittle_type = kwargs.get('brittle', None)
        peierls_type = kwargs.get('peierls', None)
        self.diff_type = diff_type
        self.disl_type = disl_type
        self.brittle_type = brittle_type
        self.peierls_type = peierls_type
        if diff_type != None and disl_type != None:
            assert(diff_type == disl_type)  # doesn't make sense to have inconsistent flow laws
        if diff_type != None:
            diffusion_creep, _ = GetRheology(diff_type)
            self.diff = diffusion_creep
        if disl_type != None:
            _, dislocation_creep = GetRheology(disl_type)
            self.disl = dislocation_creep
        if brittle_type != None:
            self.brittle = self.RheologyPrm.get_rheology(brittle_type, 'brittle')
        if self.peierls_type != None:
            self.peierls = GetPeierlsRheology(peierls_type) 
    
    
    def ReadProfile(self, file_path):
        '''
        Read a depth, pressure and temperature profile from a depth_average output
        These values are saved as class variables
        '''
        self.depths, self.pressures, self.temperatures = ReadAspectProfile(file_path)


    # todo_HK03
    def VaryWithStress(self, P, T, d, Coh, stress_range, **kwargs):
        '''
        With a set of variables, compare to the data points reported
        in experimental publications.
        Inputs:
            stress_range - a range of stress for the strain rate - stress plot
            strain_rate_range - a range of strain_rate, only affects plot
            P,d,Coh - samples of variables for the strain rate - stress plot
            kwargs:
                ax - an axis for plot
                label - label for the curve
                color - the color used for plotting
                
        '''
        assert(self.diff is not None and self.disl is not None)
        stress_range = kwargs.get("stress_range", [10.0, 1000.0])  # Mpa
        strain_rate_range = kwargs.get("strain_rate_range", None)  # s^{-1}
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        _color = kwargs.get('color', 'b')
        use_effective_strain_rate = kwargs.get('use_effective_strain_rate', True)

        stresses = np.linspace(stress_range[0], stress_range[1], 1000)
        strain_rates_diff = CreepStrainRate(self.diff, stresses, P, T, d, Coh, use_effective_strain_rate=use_effective_strain_rate)
        strain_rates_disl = CreepStrainRate(self.disl, stresses, P, T, d, Coh, use_effective_strain_rate=use_effective_strain_rate)
        strain_rates = strain_rates_diff + strain_rates_disl  # apply the isostress model
        if ax is not None:
            ax.loglog(stresses, strain_rates, '-', color=_color, label=(label + "comp"))
            ax.loglog(stresses, strain_rates_diff, '--', color=_color, label=(label + "diff"))
            ax.loglog(stresses, strain_rates_disl, '-.', color=_color, label=(label + "disl"))
        ax.set_xlabel("Stress (Mpa)")
        ax.set_xlim(stress_range[0], stress_range[1])
        if strain_rate_range is not None:
            ax.set_ylim(strain_rate_range[0], strain_rate_range[1])
        ax.set_ylabel("Strain Rate (s^-1)")
        _title = "P = %.4e, T = %.4e, d = %.4e, Coh = %.4e" % (P/1e6, T, d, Coh)
        ax.set_title(_title)

    
    def VaryWithT(self, P, stress, d, Coh, T_range, **kwargs):
        '''
        With a set of variables, compare to the data points reported
        in experimental publications.
        Inputs:
            P - a pressure to compute the strain rate
            stress - a sample of stress for the strain rate - 10^4/T plot
            d - a sample of grain size for the strain rate - 10^4/T plot
            Coh - a sample of Coh for the strain rate - 10^4/T plot
            T_range - a range of T for the strain rate - 10^4/T plot
            kwargs:
                ax - an axis for plot
                label - label for the curve
                color - the color used for plotting
        '''
        assert(self.diff is not None and self.disl is not None)
        strain_rate_range = kwargs.get("strain_rate_range", None)  # s^{-1}
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        _color = kwargs.get('color', 'b')
        use_effective_strain_rate = kwargs.get('use_effective_strain_rate', True)
        
        Ts = np.linspace(T_range[0], T_range[1], 1000)
        strain_rates_diff = CreepStrainRate(self.diff, stress, P, Ts, d, Coh, use_effective_strain_rate=use_effective_strain_rate)
        strain_rates_disl = CreepStrainRate(self.disl, stress, P, Ts, d, Coh, use_effective_strain_rate=use_effective_strain_rate)
        strain_rates = strain_rates_diff + strain_rates_disl  # apply the isostress model
        if ax is not None:
            ax.semilogy(1e4/Ts, strain_rates, '-', color=_color, label=(label + "comp"))
            ax.semilogy(1e4/Ts, strain_rates_diff, '--', color=_color, label=(label + "diff"))
            ax.semilogy(1e4/Ts, strain_rates_disl, '-.', color=_color, label=(label + "disl"))
        ax.set_xlabel("10^4 / T (K^-1)")
        ax.set_xlim(1e4 / T_range[1], 1e4 / T_range[0])
        if strain_rate_range is not None:
            ax.set_ylim(strain_rate_range[0], strain_rate_range[1])
        ax.set_ylabel("Strain Rate (s^-1)")
        _title = "P = %.4e, Stress = %.4e, d = %.4e, Coh = %.4e" % (P/1e6, stress, d, Coh)
        ax.set_title(_title)

    def MantleRheology(self, **kwargs):
        '''
        Derive mantle rheology from an aspect profile
        In this version, I would use the F factor (second invariant) as the default for computing the viscosity.
        Inputs:
            kwargs:
                rheology - type of rheology to use
                strain_rate - the strain rate used for viscosity estimation
                dEdiff, dVdiff, dAdiff_ratio, dAdisl_ratio, dEdisl, dVdisl - these factors
                    would apply a differences to the medium value in the flow law.
                save_profile - if the mantle profile of viscosity is saved as plot.
                save_json - if the derived mantle rheology is saved as a json file
                fig_path - if the save_profile is true, then a path of figure could be given
                    otherwise, a default path will be adapted.
        '''
        strain_rate = kwargs.get('strain_rate', 1e-15)
        use_effective_strain_rate = kwargs.get('use_effective_strain_rate', True)
        eta_diff = np.ones(self.depths.size)
        eta_disl = np.ones(self.depths.size)
        eta_disl13 = np.ones(self.depths.size)
        eta13 = np.ones(self.depths.size)
        eta = np.ones(self.depths.size) 
        # these options are for a differences from the central value
        dEdiff = float(kwargs.get('dEdiff', 0.0))  # numbers for the variation in the rheology
        dVdiff = float(kwargs.get('dVdiff', 0.0))
        dAdiff_ratio = float(kwargs.get("dAdiff_ratio", 1.0))
        dAdisl_ratio = float(kwargs.get("dAdisl_ratio", 1.0))
        dEdisl = float(kwargs.get('dEdisl', 0.0))
        dVdisl = float(kwargs.get('dVdisl', 0.0))
        rheology = kwargs.get('rheology', 'HK03_wet_mod')
        save_profile = kwargs.get('save_profile', 0)
        save_json = kwargs.get('save_json', 0)
        debug = kwargs.get('debug', False)
        fig_path = kwargs.get("fig_path", None)

        # First, read in the flow law and apply the difference to the medium value 
        diffusion_creep, dislocation_creep = GetRheology(rheology)
        diffusion_creep['A'] *= dAdiff_ratio
        diffusion_creep['E'] += dEdiff
        dislocation_creep['A'] *= dAdisl_ratio
        dislocation_creep['E'] += dEdisl
        diffusion_creep['V'] += dVdiff
        dislocation_creep['V'] += dVdisl
        self.diff = diffusion_creep  # record these with the class variables
        self.disl = dislocation_creep

        # Then, convert T, P as function. The T_func and P_func are used
        # in the following code to get the values
        T_func = interp1d(self.depths, self.temperatures, assume_sorted=True)
        P_func = interp1d(self.depths, self.pressures, assume_sorted=True)

        # okay, we are ready
        # Start by getting the rheology < 410 km
        depth_up = 410e3
        depth_low = 660e3
        mask_up = (self.depths < depth_up)
        eta_diff[mask_up] = CreepRheology(diffusion_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=use_effective_strain_rate)
        eta_disl[mask_up] = CreepRheology(dislocation_creep, strain_rate, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=use_effective_strain_rate)
        eta_disl13[mask_up] = CreepRheology(dislocation_creep, 1e-13, self.pressures[mask_up], self.temperatures[mask_up], use_effective_strain_rate=use_effective_strain_rate)
        eta[mask_up] = ComputeComposite(eta_diff[mask_up], eta_disl[mask_up])
        eta13[mask_up] = ComputeComposite(eta_diff[mask_up], eta_disl13[mask_up])


        # then, the rheology in the MTZ
        # Now there is no differences from the scenario we used in the upper mantle
        # in the future, more will be added.
        mask_mtz = (self.depths > depth_up) & (self.depths < depth_low)
        if True:
            # MTZ from olivine rheology
            eta_diff[mask_mtz] = CreepRheology(diffusion_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=use_effective_strain_rate)
            eta_disl[mask_mtz] = CreepRheology(dislocation_creep, strain_rate, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=use_effective_strain_rate)
            eta_disl13[mask_mtz] = CreepRheology(dislocation_creep, 1e-13, self.pressures[mask_mtz], self.temperatures[mask_mtz], use_effective_strain_rate=use_effective_strain_rate)
            eta[mask_mtz] = ComputeComposite(eta_diff[mask_mtz], eta_disl[mask_mtz])
            eta13[mask_mtz] = ComputeComposite(eta_diff[mask_mtz], eta_disl13[mask_mtz])
        
        # At last, the lower mantle
        # The diffusion creep is assumed to be the only activated mechanism in the lower mantle.
        mask_low = (self.depths > depth_low)
        jump_lower_mantle = kwargs.get('jump_lower_mantle', 30.0)
        # Computing V in the lower mantle.
        # For this, we need the T, P on the 660 boundary
        depth_lm = 660e3
        depth_max = self.depths[-1] - 10e3
        T660 = T_func(depth_lm)
        P660 = P_func(depth_lm)
        eta_diff660 = CreepRheology(diffusion_creep, strain_rate, P660, T660, use_effective_strain_rate=use_effective_strain_rate)
        # dislocation creep
        eta_disl660 = CreepRheology(dislocation_creep, strain_rate, P660, T660, use_effective_strain_rate=use_effective_strain_rate)
        eta660 = ComputeComposite(eta_diff660, eta_disl660)
        if debug:
            print("eta_diff660 = ", eta_diff660)
            print("eta_disl660 = ", eta_disl660)
            print("eta_comp660 = ", eta660)
        diff_lm = diffusion_creep.copy()
        diff_lm['V'] = 3e-6  # assign a value
        diff_lm['A'] = CreepComputeA(diff_lm, strain_rate, P660, T660, eta660*jump_lower_mantle, use_effective_strain_rate=use_effective_strain_rate)
        
        # dump json file 
        constrained_rheology = {'diffusion_creep': diffusion_creep, 'dislocation_creep': dislocation_creep, 'diffusion_lm': diff_lm}
        
        # convert aspect rheology
        print('diffusion_creep: ', diffusion_creep) # debug
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep, use_effective_strain_rate=use_effective_strain_rate)
        diffusion_lm_aspect = Convert2AspectInput(diff_lm, use_effective_strain_rate=use_effective_strain_rate)
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=use_effective_strain_rate)
        constrained_rheology_aspect = {'diffusion_creep': diffusion_creep_aspect, 'dislocation_creep': dislocation_creep_aspect, 'diffusion_lm': diffusion_lm_aspect}
        constrained_viscosity_profile = {'T': None, 'P': None, 'diffusion': None, 'dislocation': None,\
                                        'composite': None, 'dislocation_13': None,\
                                        'composite_13': None, 'depth': None}
        
        
        # lower mnatle rheology
        eta_diff[mask_low] = CreepRheologyInAspectViscoPlastic(diffusion_lm_aspect, strain_rate, self.pressures[mask_low], self.temperatures[mask_low])
        eta_disl[mask_low] = None  # this is just for visualization
        eta_disl13[mask_low] = None  # this is just for visualization
        eta[mask_low] = eta_diff[mask_low]  # diffusion creep is activated in lower mantle
        eta13[mask_low] = eta_diff[mask_low]  # diffusion creep is activated in lower mantle
        
        # Next, we visit some constraints for whole manlte rheology 
        # to see whether we match them
        # The haskel constraint
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
        if save_json == 1:
            json_path = os.path.join(RESULT_DIR, "mantle_profile_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e_dAdiff%.4e_dAdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl, dAdiff_ratio, dAdisl_ratio))
            json_path_aspect = os.path.join(RESULT_DIR, "mantle_profile_aspect_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e_dAdiff%.4e_dAdisl%.4e.json" % (rheology, dEdiff, dEdisl, dVdiff, dVdisl, dAdiff_ratio, dAdisl_ratio))
            with open(json_path, 'w') as fout:
                json.dump(constrained_rheology, fout)
            with open(json_path_aspect, 'w') as fout:
                json.dump(constrained_rheology_aspect, fout)
            print("New json: %s" % json_path)
            print("New json: %s" % json_path_aspect)
            self.output_json = json_path
            self.output_json_aspect = json_path_aspect

        # save the constrained viscosity profile
        constrained_viscosity_profile['depth'] = self.depths.copy() 
        constrained_viscosity_profile['T'] = self.temperatures.copy() 
        constrained_viscosity_profile['P'] = self.pressures.copy() 
        constrained_viscosity_profile['diffusion'] = eta_diff.copy()
        constrained_viscosity_profile['dislocation'] = eta_disl.copy()
        constrained_viscosity_profile['composite'] = eta.copy()
        constrained_viscosity_profile['composite_13'] = eta13.copy()
        constrained_viscosity_profile['dislocation_13'] = eta_disl13.copy()
        
        # plot the profile of viscosity if it's required
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
            axs[1].set_title('%s_lowerV_%.4e_haskell%.2f' % (rheology, diffusion_lm_aspect['V'], average_log_eta))
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
            if fig_path == None:
                fig_path = os.path.join(RESULT_DIR,\
                    "mantle_profile_v1_%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e_dAdiff%.4e_dAdisl%.4e.png"\
                    % (rheology, dEdiff, dEdisl, dVdiff, dVdisl, dAdiff_ratio, dAdisl_ratio))
            fig.savefig(fig_path)
            print("New figure: %s" % fig_path)
            plt.close()
            self.output_profile = fig_path
            pass
        return constrained_rheology_aspect, constrained_viscosity_profile
    
    def ConstrainRheology(self, **kwargs):
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


class PIEZOMETER():
    '''
    A class of piezometers
    Units to use in the Piezometers:
        grain size in mu m; sigma in MPa
    Attributes
    '''
    def __init__(self):
        pass

    def MehlHirth08GabbroMylonite(self, d):
        '''
        Piezometers in mylonite gabbro from the Melh Hirth 08 paper
        Inputs:
            d - grain size, mu m
        returns:
            sigma - stress, MPa
        '''
        d_min = 10.0
        d_max = 3000
        if d < d_min:
            raise ValueError("value of d smaller than the minimum limit")
        elif d > d_max:
            raise ValueError("value of d beyond the maximum limit")
        else:
            sigma = 237.0 * d**(-0.48)
        return sigma
    
    def MehlHirth08GabbroMyloniteInvert(self, sigma):
        '''
        Piezometers in mylonite gabbro from the Melh Hirth 08 paper
        Inverting the relationship, if the values of stress is beyond
        the max / min level, reture the min / max limit of the grain
        size.
        Inputs:
            sigma - stress, MPa
        returns:
            d - grain size, mu m
        '''
        sigma_min = 5.07843787
        # sigma_max = 78.4780757914
        # sigma_max = 50.0
        sigma_max = 26.0
        if sigma < sigma_min:
            d = 3000.0
        elif sigma > sigma_max:
            # d = 10.0
            # d = 35.0
            d = 100.0
        else:
            d = (sigma / 237.0) ** (-2.08333333333)
        return d


def GetRheology(rheology, **kwargs):
    '''
    read rheology parameters, and account for effects of water if it is a wet rheology
    Inputs:
        kwargs:
            dEdiff - a difference between the activation energy and the medium value in experiment
                (dVdiff, dEdisl, dVdisl) are defined in the same way
            dAdiff_ratio - a ratio of (A / A_medium) for the prefactor of the diffusion creep
                dAdisl_ratio is defined in the same way.
            use_coh - whether use the Coh or Fh2O as input into the wet rheology
    '''
    # these options are for a differences from the central value
    dEdiff = kwargs.get('dEdiff', 0.0)  # numbers for the variation in the rheology
    dVdiff = kwargs.get('dVdiff', 0.0)
    dAdiff_ratio = kwargs.get("dAdiff_ratio", 1.0)
    dAdisl_ratio = kwargs.get("dAdisl_ratio", 1.0)
    dEdisl = kwargs.get('dEdisl', 0.0)
    dVdisl = kwargs.get('dVdisl', 0.0)
    use_coh = kwargs.get("use_coh", True)
    # initiate the class object
    RheologyPrm = RHEOLOGY_PRM()
    # if the diffusion creep flow law is specified, then include it here
    if hasattr(RheologyPrm, rheology + "_diff"):
        diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
        ### if the rheology is formulated with the fugacity, convert it to using the Coh
        if diffusion_creep is not None:
            try:
                _ = diffusion_creep['wet']
            except KeyError:
                pass
            else:
                if use_coh:
                    water_creep = getattr(RheologyPrm, "water")
                    diffusion_creep['A'] = diffusion_creep['A'] / (water_creep['A'] ** diffusion_creep['r'])
                    diffusion_creep['V'] = diffusion_creep['V'] - water_creep['V'] * diffusion_creep['r']
                    diffusion_creep['E'] = diffusion_creep['E'] - water_creep['E'] * diffusion_creep['r']
            # apply the differences to the medium value
            diffusion_creep['A'] *= dAdiff_ratio
            diffusion_creep['E'] += dEdiff
            diffusion_creep['V'] += dVdiff
    else:
        diffusion_creep = None
    # if the dislocation creep flow law is specified, then include it here
    if hasattr(RheologyPrm, rheology + "_disl"):
        dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
        ### if the rheology is formulated with the fugacity, convert it to using the Coh
        if dislocation_creep is not None:
            try:
                _ = dislocation_creep['wet']
            except KeyError:
                pass
            else:
                if use_coh:
                    water_creep = getattr(RheologyPrm, "water")
                    dislocation_creep['A'] = dislocation_creep['A'] / (water_creep['A'] ** dislocation_creep['r'])
                    dislocation_creep['V'] = dislocation_creep['V'] - water_creep['V'] * dislocation_creep['r']
                    dislocation_creep['E'] = dislocation_creep['E'] - water_creep['E'] * dislocation_creep['r']
            # apply the differences to the medium value
            dislocation_creep['A'] *= dAdisl_ratio
            dislocation_creep['E'] += dEdisl
            dislocation_creep['V'] += dVdisl
    else:
        dislocation_creep = None
    # return the rheology 
    return diffusion_creep, dislocation_creep


def GetPeierlsRheology(rheology):
    '''
    read the peierls rheology parameters
    Inputs:
        rheology: a string of the type of rheology to use.
    Returns:
        peierls_creep: a dict of the flow law variables for the peierls creep
    '''
    RheologyPrm = RHEOLOGY_PRM()
    Utilities.my_assert(hasattr(RheologyPrm, rheology + "_peierls"), ValueError,\
    "The %s is not a valid option for the peierls rheology" % rheology)
    peierls_creep = getattr(RheologyPrm, rheology + "_peierls")
    return peierls_creep



def Config(_kwargs, _name, _default):
    """    def ReadProfile(self, file_path):
        self.depths, self.pressures, self.temperatures = ReadAspectProfile(file_path)ble value and assign default if not found
    """
    try:
        value = _kwargs[_name]
    except KeyError:
        value = _default
    return value


def CreepStress(creep, strain_rate, P, T, d, Coh, **kwargs):
    """
    def DislocationCreep(strain_rate, P, T, d, Coh)

    Calculate stress by flow law in form of (strain_rate / B)^(1.0 / n) * exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Mpa
    kwargs:
        use_effective_strain_rate - use the second invariant as input
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # compute F
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 2**(1/n)/3**((n+1)/2/n)
    else:
        F = 1.0
    # calculate B
    B = A * d**(-p) * Coh**r
    return F * (strain_rate / B)**(1.0 / n) * np.exp((E + P * V) / (n * R * T))


def CreepStrainRate(creep, stress, P, T, d, Coh, **kwargs):
    """
    Calculate strain rate by flow law in form of 
        B * sigma^n * exp( - (E + P * V) / (R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - stress: MPa
     - Coh: H / 10^6 Si
     - Return value: s^-1
    kwargs:
        use_effective_strain_rate - use the second invariant as input
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    r = creep['r']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    # calculate B
    # compute F
    use_effective_strain_rate = kwargs.get('use_effective_strain_rate', False)
    if use_effective_strain_rate:
        F = 3**((n+1)/2) / 2.0
    else:
        F = 1.0
    B = A * d**(-p) * Coh**r
    return F * B *stress**n * np.exp(-(E + P * V) / (R * T))


def CreepRheology(creep, strain_rate, P, T, d=1e4, Coh=1e3, **kwargs):
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
    # calculate B
    B = A * d**(-p) * Coh**r
    eta = 1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6

    return eta


def CreepComputeV(creep, strain_rate, P, T, eta, d=1e4, Coh=1e3, **kwargs):
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
        F = 1 / (2**((n-1)/n)*3**((n+1)/2/n))
    else:
        F = 1.0
    # calculate B
    B = A * d**(-p) * Coh**r
    exponential = eta / (1/2.0 * F * (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp(E / (n * R * T)) * 1e6)
    V = n * R * T * np.log(exponential) / P
    return V


def CreepComputeA(creep, strain_rate, P, T, eta, d=1e4, Coh=1e3, **kwargs):
    """
    Compute the prefactor in the rheology with other variables in a flow law (p, r, n, E, V).
    The viscosity is computed at condition of P, T and is constrained to be eta.
    Calculate viscosity by flow law in form of 0.5*(strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - creep: flow law that contains p, r, n, E, V
     - strain_rate: the strain rate to compute viscosity with
     - P: The pressure to compute viscosity, unit is Pa
     - T: The temperature to compute viscosity, unit is K
     - d: The grain size to compute viscosity, unit is mu m
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
    # calculate B
    B = (0.5*F/eta)**n * strain_rate**(1-n) * np.exp((E+P*V)/(R*T)) * (1e6)**n
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
    In this version, I am trying to take care of the F factor correctly
    Inputs:
        kwargs:
            d: um, the grain size to use, default is 1e4
            Coh: H / 10^6 Si, default is 1000.0
    Original Units in creep:
     - P: Pa
     - T: K
    Converted units in aspect_creep:
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
    d = kwargs.get('d', 1e4)
    Coh = kwargs.get('Coh', 1000.0)
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


def ComputeComposite(*Args):
    '''
    compute value of composite viscosity from value of diffusion creep and 
    dislocation creep. This will check that at least one entry is not None.
    If one of them is none, then the other entry will be directly returned
    '''
    i = 0
    indexes = []
    for Arg in Args:
        if Arg is not None:
            indexes.append(i)
        i += 1
    assert(len(indexes) > 0)  # check their is valid inputs
    if len(indexes) == 1:
        # if there is only 1 entry, just return it
        return Args[indexes[0]]
    else:
        reciprocal = 0.0
        for index in indexes:
            reciprocal += 1.0 / Args[index]
        eta_comp = 1.0 / reciprocal
        return eta_comp


#### functions for the peierls rheology

# todo_peierls
def PeierlsCreepStrainRate(creep, stress, P, T):
    """
    Calculate strain rate by flow law in form of 
        Ap * sigma^n * exp( - (E) / (R * T) * (1 - (sigma / sigmap)^p)^q)
    Units:
     - P: Pa
     - T: K
     - stress: MPa
     - Return value: s^-1
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep['A']
    p = creep['p']
    q = creep['q']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    sigp0 = creep['sigp0']
    # calculate B
    # compute F
    exponential = -(E + P*V) / (R*T) * (1 - (stress/sigp0)**p)**q
    print('exponential: ', exponential) # debug
    expo = np.exp(exponential)
    strain_rate = A * expo * stress ** n
    return strain_rate


# todo_peierls
def PeierlsCreepStress(creep, strain_rate, P, T, **kwargs):
    """
    Calculate stress by inverting the flow law in form of 
        Ap * sigma^n * exp( - (E) / (R * T) * (1 - (sigma / sigmap)^p)^q)
    Units:
     - P: Pa
     - T: K
     - strain_rate : s^-1
     - Return value: MPa
     - kwargs:
        - tolerance: the tolerance on the difference between iterations
        - iteration: number of the maximum iteration
    Pay attention to pass in the right value, this custom is inherited,
    note that the difference in the iteration is defined as the log value
    of the stress.
    """
    A = creep['A']
    p = creep['p']
    q = creep['q']
    n = creep['n']
    E = creep['E']
    V = creep['V']
    sigp0 = creep['sigp0']
    tolerance = kwargs.get('tolerance', 0.05)
    maximum_iteration = kwargs.get('iteration', 1000)
    # initialization
    difference = 1e6  # a big initial value
    stress_l = 1e-5
    stress_u = 1e12
    is_first = True
    n = 0
    while (abs(difference) > tolerance and n < maximum_iteration):
        if is_first:
            is_first = False
        else:
            # update the value of stress
            if difference > 0.0:
                stress_u = stress
            else:
                stress_l = stress
        exponential = (np.log(stress_u) + np.log(stress_l)) / 2.0
        stress = np.exp(exponential)
        strain_rate_1 = PeierlsCreepStrainRate(creep, stress, P, T)
        difference = np.log(strain_rate_1 / strain_rate)
        n += 1
    if n == maximum_iteration:
        raise(ValueError,\
            'tolerance (%f) is not reached at the end of iteration, the remanant difference is %f'\
            % (tolerance, difference))
    return stress

# todo_peierls
def PeierlsCreepRheology(creep, strain_rate, P, T, **kwargs):
    """
    Calculate stress by inverting the flow law in form of 
        Ap * sigma^n * exp( - (E) / (R * T) * (1 - (sigma / sigmap)^p)^q)
    Units:
     - P: Pa
     - T: K
     - strain_rate : s^-1
     - Return value: Pa * s
     - kwargs:
        - tolerance: the tolerance on the difference between iterations
        - iteration: number of the maximum iteration
    Pay attention to pass in the right value, this custom is inherited,
    note that the difference in the iteration is defined as the log value
    of the stress.
    """
    tolerance = kwargs.get('tolerance', 0.05)
    maximum_iteration = kwargs.get('iteration', 1000)
    stress = PeierlsCreepStress(creep, strain_rate, P, T, iteration=maximum_iteration, tolerance=tolerance)
    eta = 1e6 * stress / 2.0 / strain_rate
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
    Operator.ConstrainRheology(save_profile=save_profile, include_lower_mantle=include_lower_mantle)


# todo_HK03
def DeriveMantleRheology(file_path, **kwargs):
    '''
    Derive a Mantle rheology profile following certain procedures
    Inputs:
        file_path(str): a profile from ASPECT
    '''
    mantle_rheology_scheme = kwargs.get("rheology", "HK03_wet_mod")
    use_effective_strain_rate = kwargs.get("use_effective_strain_rate", True)
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
    print("use_effective_strain_rate: ", use_effective_strain_rate)
    fig_dir = os.path.join(RESULT_DIR, 
                    "%s_dEdiff%.4e_dEdisl%.4e_dVdiff%4e_dVdisl%.4e_dAdiff%.4e_dAdisl%.4e"\
                    % (mantle_rheology_scheme, dEdiff, dEdisl, dVdiff, dVdisl, dAdiff_ratio, dAdisl_ratio))
    if os.path.isdir(fig_dir):
        rmtree(fig_dir)
    os.mkdir(fig_dir)
    fig_path = os.path.join(fig_dir, "mantle_profile.png")
    if version == 0:
        Operator.MantleRheology(rheology=mantle_rheology_scheme,save_profile=save_profile,\
        use_effective_strain_rate=use_effective_strain_rate, save_json=True, fig_path=fig_path)
    elif version == 1:
        Operator.MantleRheology(rheology=mantle_rheology_scheme,save_profile=save_profile,\
            dEdiff=dEdiff, dEdisl=dEdisl, dVdiff=dVdiff, dVdisl=dVdisl,\
                dAdiff_ratio=dAdiff_ratio, dAdisl_ratio=dAdisl_ratio,
                use_effective_strain_rate=use_effective_strain_rate, save_json=True, fig_path=fig_path)
    else:
        raise CheckValueError('%d is not a valid version of Mantle Rheology' % version)

    fig = plt.figure(tight_layout=True, figsize=(8, 10))
    
    # strain rate vs stress
    gs = gridspec.GridSpec(2, 1)
    ax = fig.add_subplot(gs[0, 0])
    # plot reference point
    # note that the option of use_effective_strain_rate
    # is set as False to compare with the HK03 experiments
    PlotHK03DataFig2(ax, "wet")
    P = 300e6 # Pa
    T = 1250 + 273.15 # K
    d = 15 # um
    coh = 125.0 # / 10^ Si
    stress_range = [10.0, 1000.0]  # MPa
    strain_rate_range = [1e-7, 1e-1]  # s^-1, for plotting
    Operator.VaryWithStress(P, T, d, coh, stress_range, ax=ax,\
                            strain_rate_range=strain_rate_range, use_effective_strain_rate=False,\
                            label="")
    ax.legend()
    # strain rate vs 10^4/T
    ax = fig.add_subplot(gs[1, 0])
    # plot reference point, from figure 3
    PlotHK03DataFig3(ax, "wet") 
    P = 100e6 # Pa
    stress = 50.0 # Mpa
    d = 15 # um
    coh = 40.0 # / 10^ Si
    strain_rate_range = [1e-7, 1e-3]  # s^-1, for plotting
    T_range = [1000.0 + 273.15, 1400.0 + 273.15]
    Operator.VaryWithT(P, stress, d, coh, T_range, ax=ax,\
                            strain_rate_range=strain_rate_range, use_effective_strain_rate=False,
                            label="")
    ax.legend()
    fig_path = os.path.join(fig_dir, 'fit_HK03_fig2.png')
    fig.savefig(fig_path)
    print("New figure: ", fig_path)  # screen output


def PlotHK03DataFig2(ax, _type, **kwargs):
    '''
    Inputs:
        ax - an axis to plot
        _type - "dry" or "wet"
        kwargs:
            color - color of plotting
    '''
    if _type == 'dry':
        raise NotImplementedError()
    if _type == 'wet':
        file_path = os.path.join(ASPECT_LAB_DIR, "files", "ref_data", "HK03_fig2b_comp")
        file_path_disl = os.path.join(ASPECT_LAB_DIR, "files", "ref_data", "HK03_fig2b_disl")
    assert(os.path.isfile(file_path))
    assert(os.path.isfile(file_path_disl))
    _color = kwargs.get("color", 'b')

    # composite strain rate
    data = np.loadtxt(file_path)
    stress = 10**data[:, 0]  # MPa
    strain_rate = 10**data[:, 1] # s^-1
    ax.loglog(stress, strain_rate, '*', color=_color)
    # dislocation strain rate
    data = np.loadtxt(file_path_disl)
    stress = 10**data[:, 0]  # MPa
    strain_rate = 10**data[:, 1] # s^-1
    ax.loglog(stress, strain_rate, 'o', color=_color)


def PlotHK03DataFig3(ax, _type, **kwargs):
    '''
    Inputs:
        ax - an axis to plot
        _type - "dry" or "wet"
        kwargs:
            color - color of plotting
    '''
    if _type == 'dry':
        raise NotImplementedError()
    if _type == 'wet':
        file_path_diff = os.path.join(ASPECT_LAB_DIR, "files", "ref_data", "HK03_fig3b_diff_100")
        file_path_disl = os.path.join(ASPECT_LAB_DIR, "files", "ref_data", "HK03_fig3b_disl_100")
    assert(os.path.isfile(file_path_diff))
    assert(os.path.isfile(file_path_disl))
    _color = kwargs.get("color", 'b')

    # diffusion strain rate
    data = np.loadtxt(file_path_diff)
    Ts = 1e4 / data[:, 0]  # K
    strain_rate = 10**data[:, 1] # s^-1
    ax.semilogy(1e4 / Ts, strain_rate, 's', color=_color, label="diffusion, P = 100 MPa")
    # dislocation strain rate
    data = np.loadtxt(file_path_disl)
    Ts = 1e4 / data[:, 0]  # K
    strain_rate = 10**data[:, 1] # s^-1
    ax.semilogy(1e4 / Ts, strain_rate, 'o', color=_color, label="dislocation, P = 100 MPa")

###
# functions for deriving the strenght profile
###
class STRENGTH_PROFILE(RHEOLOGY_OPR):

    def __init__(self, **kwargs):
        RHEOLOGY_OPR.__init__(self)
        self.Sigs = None
        self.Sigs_brittle = None
        self.etas_brittle = None
        self.Sigs_viscous = None
        self.etas_viscous = None
        self.etas_peierls = None
        self.Zs = None
        self.Etas = None
        self.Computed = False
        # todo_peierls
        self.peierls_type = None
        self.max_depth = kwargs.get('max_depth', 80e3)
        self.T_type = kwargs.get("T_type", "hpc")
        self.SetRheologyByName(diff=kwargs.get("diff", None),\
                               disl=kwargs.get("disl", None),\
                               brittle=kwargs.get("brittle", None),\
                               peierls=kwargs.get("peierls", None))

    def Execute(self, **kwargs):
        '''
        Compute the strength profile
        '''
        year = 365 * 24 * 3600.0
        compute_second_invariant = kwargs.get('compute_second_invariant', False)
        brittle = self.brittle
        assert(self.diff_type != None or self.disl_type != None)
        assert(brittle != None)
        averaging = kwargs.get('averaging', 'harmonic')
        strain_rate = kwargs.get('strain_rate', 1e-14)
        # rheology_prm = RHEOLOGY_PRM()
        # self.plastic = rheology_prm.ARCAY17_plastic
        # self.dislocation_creep = rheology_prm.ARCAY17_disl
        Zs = np.linspace(0.0, self.max_depth, 100)
        if self.T_type == 'hpc':
            # use a half space cooling
            Ts = temperature_halfspace(Zs, 40e6*year, Tm=1573.0) # adiabatic temperature
        elif self.T_type == 'ARCAY17':
            Ts = 713 * Zs / 78.245e3  + 273.14# geotherm from Arcay 2017 pepi, figure 3d 2
        else:
            raise NotImplementedError
        Tliths = temperature_halfspace(Zs, 40e6*year, Tm=1573.0) # adiabatic temperature
        Ps = pressure_from_lithostatic(Zs, Tliths)
        # brittle self.brittle
        if brittle["type"] == "stress dependent":
            # note this is questionable, is this second order invariant
            self.Sigs_brittle = StressDependentYielding(Ps, brittle["cohesion"], brittle["friction"], brittle["ref strain rate"], brittle["n"], strain_rate)
        elif brittle["type"] == "Coulumb":
            self.Sigs_brittle = CoulumbYielding(Ps, brittle["cohesion"], brittle["friction"])
        elif brittle["type"] == "Byerlee":
            self.Sigs_brittle = Byerlee(Ps)
        else:
            raise NotImplementedError()
        self.etas_brittle = self.Sigs_brittle / 2.0 / strain_rate
        # viscous stress
        # Note on the d and coh:
        #      disl - not dependent on d;
        #      coh - the one in the Arcay paper doesn't depend on Coh
        etas_diff = None
        etas_disl = None
        etas_peierls = None
        if self.diff_type is not None:
            etas_diff = CreepRheology(self.diff, strain_rate, Ps, Ts,\
                                          1e4, 1000.0, use_effective_strain_rate=compute_second_invariant)
        if self.disl_type is not None:
            etas_disl = CreepRheology(self.disl, strain_rate, Ps, Ts,\
                                          1e4, 1000.0, use_effective_strain_rate=compute_second_invariant)
        if self.peierls_type is not None:
            self.etas_peierls = np.zeros(Ts.size)
            for i in range(Ts.size):
                P = Ps[i]
                T = Ts[i]
                self.etas_peierls[i] = PeierlsCreepRheology(self.peierls, strain_rate, P, T)
            self.Sigs_peierls = 2.0 * strain_rate * self.etas_peierls
        self.etas_viscous = ComputeComposite(etas_diff, etas_disl)
        self.Sigs_viscous = 2.0 * strain_rate * self.etas_viscous
        # Sigs_viscous = CreepStress(self.creep, strain_rate, Ps, Ts, 1e4, 1000.0) # change to UI
        if averaging == 'harmonic':
            self.Etas = ComputeComposite(self.etas_viscous, self.etas_brittle, self.etas_peierls)
        else:
            raise NotImplementedError()
        self.Sigs = 2 * strain_rate * self.Etas
        self.Zs = Zs
        self.computed = True

    def PlotStress(self, **kwargs):
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        label_components = kwargs.get('label_components', False)
        plot_stress_by_log = kwargs.get('plot_stress_by_log', False)
        if label_components:
            label_brittle = "brittle"
            label_viscous = "viscous"
            label_peierls = "peierls"
        else:
            label_brittle = None
            label_viscous = None
            label_peierls = None
        _color = kwargs.get('color', 'b')
        if ax == None:
            raise NotImplementedError()
        # make plots
        mask = (self.Etas > 1e-32) # get the reasonable values, the peierls creep may return inf values
        if plot_stress_by_log:
            # plot the log values of the stress on the x axis
            ax.semilogx(self.Sigs[mask]/1e6, self.Zs[mask]/1e3, color=_color, label=label)
            ax.semilogx(self.Sigs_brittle[mask]/1e6, self.Zs[mask]/1e3, '.', color=_color, label=label_brittle)
            ax.semilogx(self.Sigs_viscous[mask]/1e6, self.Zs[mask]/1e3, '--', color=_color, label=label_viscous)
            if self.peierls_type is not None:
                ax.semilogx(self.Sigs_peierls[mask]/1e6, self.Zs[mask]/1e3, '-.', color=_color, label=label_peierls)
            ax.set_xlim([0.0, 10**(3.5)])
        else:
            # plot the log values of the stress on the x axis
            ax.plot(self.Sigs[mask]/1e6, self.Zs[mask]/1e3, color=_color, label=label)
            ax.plot(self.Sigs_brittle[mask]/1e6, self.Zs[mask]/1e3, '.', color=_color, label=label_brittle)
            ax.plot(self.Sigs_viscous[mask]/1e6, self.Zs[mask]/1e3, '--', color=_color, label=label_viscous)
            if self.peierls_type is not None:
                ax.plot(self.Sigs_peierls[mask]/1e6, self.Zs[mask]/1e3, '-.', color=_color, label=label_peierls)
            x_max = np.ceil(np.max(self.Sigs[mask]/1e6) / 100.0) * 100.0
            ax.set_xlim([0.0, x_max])
        ax.set_xlabel("Second invariant of the stress tensor (MPa)")
        ax.set_ylabel("Depth (km)")
    
    def PlotViscosity(self, **kwargs):
        ax = kwargs.get('ax', None)
        label = kwargs.get('label', None)
        label_components = kwargs.get('label_components', False)
        if label_components:
            label_brittle = "brittle"
            label_viscous = "viscous"
            label_peierls = "peierls"
        else:
            label_brittle = None
            label_viscous = None
            label_peierls = None
        _color = kwargs.get('color', 'b')
        if ax == None:
            raise NotImplementedError()
        # plot viscosity
        mask = (self.Etas > 1e-32) # get the reasonable values, the peierls creep may return inf values
        ax.semilogx(self.Etas[mask], self.Zs[mask]/1e3, color=_color, label=label)
        ax.semilogx(self.etas_brittle[mask], self.Zs[mask]/1e3, '.', color=_color, label=label_brittle)
        ax.semilogx(self.etas_viscous[mask], self.Zs[mask]/1e3, '--', color=_color, label=label_viscous)
        if self.peierls_type is not None:
            ax.semilogx(self.etas_peierls[mask], self.Zs[mask]/1e3, '-.', color=_color, label=label_peierls)
        ax.set_xlabel("Viscosity (Pa * s)")
        ax.set_ylabel("Depth (km)")
        x_min = 1e18
        x_max = 1e24
        ax.set_xlim([x_min, x_max])
        

def PlotStrengh(Operator, fig_path_base, **kwargs):
    '''
    Plot the shear zone strenght profile
    kwargs:
        compute_second_invariant: whether using the second invariant as the stress.
    '''
    compute_second_invariant = kwargs.get('compute_second_invariant', False)
    strain_rate = kwargs.get('strain_rate', None)
    if strain_rate is None:
        strain_rates = [1e-13, 1e-14, 1e-15]
    else:
        strain_rates = [strain_rate]
    fig = plt.figure(tight_layout=True, figsize=[10, 10])
    gs = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    colors = ['b', 'g', 'r']
    # 1e-13 
    i = 0
    for strain_rate in strain_rates:
        _color = colors[i]
        Operator.Execute(creep_type='disl', strain_rate=strain_rate,\
        compute_second_invariant=compute_second_invariant)
        # plot stress
        label = "Strain Rate = %.1e" % strain_rate
        Operator.PlotStress(ax=ax0, color=_color, label=label, label_components=(i==0))
        Operator.PlotStress(ax=ax1, color=_color, plot_stress_by_log=True)
        # plot viscosity
        Operator.PlotViscosity(ax=ax2, color=_color)
        i += 1
    ax0.invert_yaxis()
    ax0.legend()
    ax1.invert_yaxis()
    ax1.legend()
    ax2.invert_yaxis()
    ax2.legend()
    # title, add the name of the rheology being used
    fig_title = "brittle type: " + Operator.brittle_type
    if Operator.diff_type is not None:
        fig_title += ', creep (diff): ' + Operator.diff_type
    if Operator.disl_type is not None:
        fig_title += ', creep (disl): ' + Operator.disl_type
    if Operator.peierls_type is not None:
        fig_title += ', peierls: ' + Operator.peierls_type
    fig.suptitle(fig_title)
    # figure path
    fig_path = fig_path_base.split('.')[0]
    fig_path += ('_brittle_' + Operator.brittle_type)
    try:
        if Operator.diff_type is not None:
            fig_path += ("_diff_" + Operator.diff_type)
        if Operator.disl_type is not None:
            fig_path += ("_disl_" + Operator.disl_type)
        if Operator.peierls_type is not None:
            fig_path += ("_peierls_" + Operator.peierls_type)
    except TypeError:
        # there is no "disl_type" for the operator
        fig_path = fig_path_base
    fig_path += ("_T_" + Operator.T_type)
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


# todo_r_json
def PlotStrainRateStress(diff, disl, dA_diff_ratio, dE_diff, dV_diff,\
                        dA_disl_ratio, dE_disl, dV_disl, grain_size, coh,\
                        fh2o, use_coh, **kwargs):
    '''
    Plot a strain rate - stress profile
    using a json file as input.
    Inputs:
        dE_diff - a difference between the activation energy and the medium value in experiment
            (dV_diff, dE_disl, dV_disl) are defined in the same way
        dA_diff_ratio - a ratio of (A / A_medium) for the prefactor of the diffusion creep
            dA_disl_ratio is defined in the same way.
        kwargs:
            color - the color to plot with
    '''
    ax = kwargs['ax']
    # color of the plot
    _color = kwargs.get('color', "tab:blue")
    # options for the range of stress
    stress_min = 5.0 # MPa
    stress_max = 500.0 # MPa
    stress_N = 100  # integer, number of sample points for the stress
    # options for the range of strain rate 
    strain_rate_min = 10**(-15)  # s^-1
    strain_rate_max = 10**(-10) # s^-1
    # options for the temperature
    T = kwargs.get("T", 875.0)  # C
    # options for the pressure 
    P = 2500 * 10 * 3000.0 # Pa, has to be the same shape as Ts
    # options of rheology 
    # include the differences to the reference valuea
    wet = False

    if use_coh:
        wet_variable = coh
    else:
        wet_variable = fh2o
    if diff != "":
        rheology_diff, _  = GetRheology(diff,
            dEdiff=dE_diff, dVdiff=dV_diff, dAdiff_ratio=dA_diff_ratio,
            use_coh=use_coh)
        if rheology_diff['r'] > 1e-6:
            wet = True
    if disl != "":
        _, rheology_disl = GetRheology(disl,
            dAdisl_ratio=dA_disl_ratio, dEdisl=dE_disl, dVdisl=dV_disl,
            use_coh=use_coh)
        if rheology_disl['r'] > 1e-6:
            wet = True
    # initiation 
    stresses = 10.0**(np.linspace(np.log10(stress_min), np.log10(stress_max), stress_N))
    # plot
    # compute strain rates
    strain_rates = np.zeros(stresses.shape)
    for j in range(stress_N):
        stress = stresses[j]
        # for olivine
        strain_rate = 0.0
        if diff != "":
            strain_rate_diff = CreepStrainRate(rheology_diff, stress, P, T + 273.15, grain_size, wet_variable)  # dry rheology, Coh is not dependent
            strain_rate += strain_rate_diff
        if disl != "":
            strain_rate_disl = CreepStrainRate(rheology_disl, stress, P, T + 273.15, grain_size, wet_variable)  # dry rheology, Coh is not dependent
            strain_rate += strain_rate_disl
        strain_rates[j] =  strain_rate # the iso stress model
    # plot stress vs strain_rate
    ax.loglog(stresses, strain_rates, color=_color)
    # add plot options
    ax.set_xlabel("Stress (MPa)")
    ax.set_ylabel('Strain Rate (s^-1)')
    ax.set_xlim([stress_min, stress_max])
    ax.set_ylim([strain_rate_min, strain_rate_max])
    ax.set_title("%.2e C, %.2e Pa" % (T, P))
    ax.grid(True)
    # return the label and patch
    label = "d: " + str(grain_size)
    # append value of coh in case of wet rheology
    if wet:
        if use_coh:
            label += ", coh: " + str(coh)
        else:
            label += ", fh2o: " + str(fh2o)
    label += '\n'
    if diff != "":
        label += "diff: " + str(diff) + ", " + "\n" + str(rheology_diff) + "\n"
    if disl != "":
        label += ", disl: " + str(disl) + ", " + "\n" + str(rheology_disl) + "\n"
    patch = mpatches.Patch(color=_color)
    r_dict = {'label': label, 'patch': patch}
    return r_dict


def PlotViscosityTemperature(diff, disl, dA_diff_ratio, dE_diff, dV_diff,\
                            dA_disl_ratio, dE_disl, dV_disl, grain_size, coh,\
                            fh2o, use_coh, **kwargs):
    '''
    Plot the viscosity vs temperature
    Inputs:
        dE_diff - a difference between the activation energy and the medium value in experiment
            (dV_diff, dE_disl, dV_disl) are defined in the same way
        dA_diff_ratio - a ratio of (A / A_medium) for the prefactor of the diffusion creep
            dA_disl_ratio is defined in the same way.
        kwargs:
            ax - the axis to plot with
            color - the color to plot with
    '''
    ax = kwargs['ax']
    strain_rate = kwargs.get('strain_rate', 1e-15)
    # color of the plot
    _color = kwargs.get('color', "tab:blue")
    # options for temperature
    T_min = 0.0 # C
    T_max = 1000 # C
    T_N = 1000
    # options for viscosity (only affects plotting)
    eta_min = 1e17
    eta_max = 1e25
    # Options for the pressure
    P = 1.6e9 # 50 km
    Ts = np.linspace(T_min, T_max, T_N)
    # options of rheology 
    # include the differences to the reference value
    wet = False
    
    if use_coh:
        wet_variable = coh
    else:
        wet_variable = fh2o
    if diff != "":
        rheology_diff, _  = GetRheology(diff,
            dEdiff=dE_diff, dVdiff=dV_diff, dAdiff_ratio=dA_diff_ratio, use_coh=use_coh)
        if rheology_diff['r'] > 1e-6:
            wet = True
    if disl != "":
        _, rheology_disl = GetRheology(disl,
            dAdisl_ratio=dA_disl_ratio, dEdisl=dE_disl, dVdisl=dV_disl, use_coh=use_coh)
        if rheology_disl['r'] > 1e-6:
            wet = True
    grain_size = grain_size
    # initiate and compute viscosity
    etas = np.zeros(T_N)
    for j in range(T_N):
        T = Ts[j] + 273.15     
        # option of the dry olivine
        eta_diff = None
        eta_disl = None
        if diff != "":
            eta_diff =  CreepRheology(rheology_diff, strain_rate, P, T, grain_size, wet_variable)
        if disl != "":
            eta_disl =  CreepRheology(rheology_disl, strain_rate, P, T, grain_size, wet_variable)
        etas[j] = ComputeComposite(eta_diff, eta_disl)
    ax.semilogy(Ts, etas, color=_color)
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel('Viscosity (Pa s)')
    ax.set_xlim([T_min, T_max])
    ax.set_ylim([eta_min, eta_max])
    ax.set_title("%.2e s^-1, %.2e Pa" % (strain_rate, P))
    ax.grid(True)
    # return the label and patch
    label = "d: " + str(grain_size)
    # append value of coh in case of wet rheology
    if wet:
        if use_coh:
            label += ", coh: " + str(coh)
        else:
            label += ", fh2o: " + str(fh2o)
    label += '\n'
    if diff != "":
        label += "diff: " + str(diff) + ", " + "\n" + str(rheology_diff) + "\n"
    if disl != "":
        label += ", disl: " + str(disl) + ", " + "\n" + str(rheology_disl) + "\n"
    patch = mpatches.Patch(color=_color)
    r_dict = {'label': label, 'patch': patch}
    return r_dict


def PlotRheologySummaryJson(json_file):
    '''
    Plot the summary of rheology
    Inputs:
        json_file: path to the json_file
    '''
    assert(os.path.isfile(json_file))
    # parse the inputs from the json file
    RheologyJson = RHEOLOGY_JSON()
    RheologyJson.read_json(json_file)
    RheologyJson.check()
    # inputs from the plot features
    plots = RheologyJson.GetPlotFeatures()
    n_plots = len(plots)
    unit_size = 5
    n_legend_plot = int(np.ceil(len(RheologyJson.GetRheologyFeatures()) / 5.0))# these are just legends
    n_plot_col = 2
    n_plot_row = int(np.ceil(n_plots / n_plot_col)) + n_legend_plot
    # options of colors for different rheologies
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray'] 
    # inputs from the rheology features
    rheologies = RheologyJson.GetRheologyFeatures()
    assert(len(rheologies) <= 8)  # maximum number of lines
    # initiate the plot
    fig = plt.figure(tight_layout=True, figsize=(n_plot_col*unit_size, n_plot_row*unit_size))
    gs = gridspec.GridSpec(n_plot_row, n_plot_col)
    # plot figures
    i_plot = 0
    labels = []
    patches = []
    for plotOpt in plots:
        # figure out the col and row
        plot_col = i_plot % n_plot_col
        plot_row = int(np.floor(i_plot / n_plot_col)) + n_legend_plot
        # get the axis & values to plot
        ax = fig.add_subplot(gs[plot_row, plot_col])
        strain_rates = plotOpt.GetStrainRates()
        plot_type = plotOpt.GetType()
        T = plotOpt.GetT()
        i = 0
        # loop for rheologies
        for rheologyOpt in rheologies:
            # get the value of grain size
            # and strain rates for the plot
            # either 1 or multiple values will
            # be used for each rheology we
            # want to plot
            if len(strain_rates) == 1:
                strain_rate = strain_rates[0] 
            else:
                strain_rate = strain_rates[i]
            _color = colors[i]
            if plot_type == "strain rate vs stress":
                r_dict = PlotStrainRateStress(*rheologyOpt.to_RheologyInputs(), ax=ax, color=_color, T=T)
            elif plot_type == "viscosity vs temperature":
                PlotViscosityTemperature(*rheologyOpt.to_RheologyInputs(), ax=ax, strain_rate=strain_rate, color=_color)
            # append the label & patch if this is the first plot
            if i_plot == 0:
                labels.append(r_dict['label'])
                patches.append(r_dict['patch'])
            i += 1
        i_plot += 1
    # add a plot of labels on top
    ax = fig.add_subplot(gs[0:n_legend_plot, :])
    ax.legend(patches, labels, loc='upper left', frameon=False, fontsize=8) 
    fig_path = RheologyJson.GetFigurePath()
    fig.savefig(fig_path)
    print("Save figure: ", fig_path)
    
    #
    json_output_path = fig_path.split('.')[0] + ".json"
    copy(json_file, json_output_path)
    print("Save json file: ", json_output_path)


def PlotShearZoneRheologySummary(**kwargs):
    '''
    A function to plot the shear zone rheology, in a strain rate vs stress plot.
    Different rocks / minerals are considered with their flow laws.
    Options of grain sizes are also included.
    An example is the Figure 13 in the Mehl and Hirth 2008 paper.
    kwargs:
        option - option of rheology to use
            none - use the value presented in experiments
            Q_503_gabbro - use a 503 kj /mol in the gabbro aggregate rheology
            Q_473_and_A_gabbro - use a 473 kj /mol and a smaller A in the gabbro aggregate rheology
            d_100_gabbro - use a 100 mu m grain size in gabbro aggregate rheology
    '''
    option = kwargs.get('option', None)
    if option is not None:
        assert(option in ['Q_503_gabbro', 'Q_473_and_A_gabbro', 'd_100_gabbro'])
    # options for rheologies
    rheology_olivine = "HK03_dry"
    rheology_wet_olivine = "HK03_dry"
    rheology_dry_anorthite = "Rybachi_2000_An100_dry"
    rheology_gabbro_mylonite = "Dimanov_Dresen_An50Di35D_dry"
    rheology_basalt = "ST1981_basalt"
    Piezometer = PIEZOMETER()  # a piezometer object
    ### First part, plot strain rate vs stress (e.g. Mehl & Hirth 2008)
    # options for stress
    stress_min = 5.0 # MPa
    stress_max = 500.0 # MPa
    stress_N = 100  # integer, number of sample points for the stress
    # options for the strain rate, only affect ploting
    strain_rate_min = 10**(-15)  # s^-1
    strain_rate_max = 10**(-10) # s^-1
    # options for the temperature, pressure
    Ts = np.array([875.0, 925.0])  # C
    T_N = Ts.size
    depths = np.array([2500, 3000])
    assert(T_N == depths.size)
    Ps = depths * 10 * 3000.0 # Pa, has to be the same shape as Ts
    assert(T_N == Ps.size)
    # options for the olivine rheology
    dry_olivine_diff, dry_olivine_disl = GetRheology(rheology_olivine)
    d_olivine_dry = 3000  # mu m
    # options for the wet olivine rheology
    # note that we use the rheology in ASPECT here (in UI)
    source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    d_olivine_wet = 10000  # mu m
    Operator = RHEOLOGY_OPR()
    Operator.ReadProfile(da_file) # read profile
    rheology_aspect, _ = Operator.MantleRheology(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=30e3,\
    dVdiff=-5.5e-6, dVdisl=2.12e-6, save_profile=1, dAdiff_ratio=0.33333247873, dAdisl_ratio=1.040297619,\
    jump_lower_mantle=15.0)
    wet_olivine_diff_aspect = rheology_aspect['diffusion_creep']
    wet_olivine_diff_aspect['d'] = d_olivine_wet / 1e6 # use m
    wet_olivine_disl_aspect = rheology_aspect['dislocation_creep']
    wet_olivine_disl_aspect['d'] = d_olivine_wet / 1e6 # use m
    # options for the dry anorthite rheology
    dry_anorthite_diff, dry_anorthite_disl = GetRheology(rheology_dry_anorthite)
    d_anorthite_mylonite_dry_array = [100.0, 300.0]
    # options for the gabbro rheology
    gabbro_diff, gabbro_disl = GetRheology(rheology_gabbro_mylonite)
    gabbro_small_p_diff, gabbro_small_p_disl = GetRheology(rheology_gabbro_mylonite)
    # modify to the original flow law
    if option is None:
        d_gabbro_mylonite_dry_array = [35.0, 100.0]
    elif option == 'Q_503_gabbro':
        # use a big value of activation energy
        gabbro_diff['E'] = 503e3  # test
        gabbro_small_p_diff['E'] = 503e3  # test
        d_gabbro_mylonite_dry_array = [35.0, 100.0]
    elif option == 'Q_473_and_A_gabbro':
        # use both a higher value of E and 
        # a small value of A
        gabbro_diff['A'] = 2.14375e12  # test
        gabbro_diff['E'] = 473e3  # test
        gabbro_small_p_diff['A'] = 2.14375e12  # test
        gabbro_small_p_diff['E'] = 473e3  # test
        d_gabbro_mylonite_dry_array = [35.0, 100.0]
    elif option == 'd_100_gabbro':
        # use big value of d, even for colder temperature
        d_gabbro_mylonite_dry_array = [100.0, 100.0]
    gabbro_small_p_diff['p'] = 2.4  # test
    # d_gabbro_mylonite_dry = 35.0  # mu m
    d_gabbro_coarse = 3000.0  # mu m
    # options for the basalt rheology
    basalt_diff, basalt_disl = GetRheology(rheology_basalt)
    d_basalt_coarse = 3000.0 # mu m
    ## initializing
    fig = plt.figure(tight_layout=True, figsize = (8 * T_N, 10))
    gs = gridspec.GridSpec(2, T_N)
    stresses = 10.0**(np.linspace(np.log10(stress_min), np.log10(stress_max), stress_N))
    for i in range(Ts.size):
        # initalize for a subplot
        T = Ts[i]
        P = Ps[i]
        depth = depths[i]
        ax = fig.add_subplot(gs[0, i])
        strain_rate_dict = {}  # use for storing all the stress results
        grain_size_dict = {}  # use for stroing all the grain size results
        strain_rate_array_olivine = np.zeros(stresses.shape)
        strain_rate_array_anorthite = np.zeros(stresses.shape)
        strain_rate_array_gabbro_mylonite = np.zeros(stresses.shape)
        strain_rate_array_gabbro_mylonite_small_p = np.zeros(stresses.shape)
        strain_rate_array_gabbro_coarse = np.zeros(stresses.shape)
        strain_rate_array_basalt = np.zeros(stresses.shape)
        for j in range(stress_N):
            stress = stresses[j]
            # for olivine
            strain_rate_diff = CreepStrainRate(dry_olivine_diff, stress, P, T + 273.15, d_olivine_dry, 0.0)  # dry rheology, Coh is not dependent
            strain_rate_disl = CreepStrainRate(dry_olivine_disl, stress, P, T + 273.15, d_olivine_dry, 0.0)  # dry rheology, Coh is not dependent
            strain_rate_array_olivine[j] =  strain_rate_diff + strain_rate_disl # the iso stress model
            # for dry anorthite
            d_anorthite_mylonite_dry = d_anorthite_mylonite_dry_array[i]  # use a constant grain size
            strain_rate_diff = CreepStrainRate(dry_anorthite_diff, stress, P, T + 273.15, d_anorthite_mylonite_dry, 0.0)  # dry rheology, Coh is not dependent
            strain_rate_disl = CreepStrainRate(dry_anorthite_disl, stress, P, T + 273.15, d_anorthite_mylonite_dry, 0.0)  # dry rheology, Coh is not dependent
            strain_rate_array_anorthite[j] =  strain_rate_diff + strain_rate_disl # the iso stress model
            # d_gabbro_mylonite_dry = Piezometer.MehlHirth08GabbroMyloniteInvert(stress)  # use piezometer for gabbro mylonite
            # gabbro and gabbro with a smaller p
            d_gabbro_mylonite_dry = d_gabbro_mylonite_dry_array[i]  # use a constant grain size
            strain_rate_diff = CreepStrainRate(gabbro_diff, stress, P, T + 273.15, d_gabbro_mylonite_dry, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_disl = CreepStrainRate(gabbro_disl, stress, P, T + 273.15, d_gabbro_mylonite_dry, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_array_gabbro_mylonite[j] =  strain_rate_diff + strain_rate_disl # the iso stress model
            # print("strain_rate_diff: ", strain_rate_diff, ", strain_rate_disl: ", strain_rate_disl)  # debug
            strain_rate_diff = CreepStrainRate(gabbro_small_p_diff, stress, P, T + 273.15, d_gabbro_mylonite_dry, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_disl = CreepStrainRate(gabbro_small_p_disl, stress, P, T + 273.15, d_gabbro_mylonite_dry, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_array_gabbro_mylonite_small_p[j] =  strain_rate_diff + strain_rate_disl # the iso stress model
            # coarse gabbro
            strain_rate_diff = CreepStrainRate(gabbro_diff, stress, P, T + 273.15, d_gabbro_coarse, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_disl = CreepStrainRate(gabbro_disl, stress, P, T + 273.15, d_gabbro_coarse, 0.0)  # dry gabbro rheology, Coh is not dependent
            strain_rate_array_gabbro_coarse[j] =  strain_rate_diff + strain_rate_disl # the iso stress model
            # for basalt, only dislocation
            strain_rate_disl = CreepStrainRate(basalt_disl, stress, P, T + 273.15, d_basalt_coarse, 0.0)  # basalt rheology
            strain_rate_array_basalt[j] = strain_rate_disl
        strain_rate_dict['olivine_dry'] = strain_rate_array_olivine
        strain_rate_dict['anorthite_mylonite_dry'] = strain_rate_array_anorthite
        strain_rate_dict['gabbro_mylonite'] = strain_rate_array_gabbro_mylonite
        strain_rate_dict['gabbro_mylonite_small_p'] = strain_rate_array_gabbro_mylonite_small_p
        strain_rate_dict['gabbro_course'] = strain_rate_array_gabbro_coarse
        strain_rate_dict['basalt'] = strain_rate_array_basalt
        grain_size_dict['olivine_dry'] = d_olivine_dry
        grain_size_dict['anorthite_mylonite_dry'] = d_anorthite_mylonite_dry
        grain_size_dict['gabbro_mylonite'] = d_gabbro_mylonite_dry
        grain_size_dict['gabbro_mylonite_small_p'] = d_gabbro_mylonite_dry
        grain_size_dict['gabbro_course'] = d_gabbro_coarse
        grain_size_dict['basalt'] = d_basalt_coarse
        for key, value in strain_rate_dict.items():
            _label = key + ", d = " + str(grain_size_dict[key]) + " um"
            ax.loglog(stresses, value, label=_label)
        ax.legend()
        ax.set_xlabel("Stress (MPa)")
        ax.set_xlim([stress_min, stress_max])
        ax.set_ylim([strain_rate_min, strain_rate_max])
        ax.set_ylabel('Strain Rate (s^-1)')
        ax.set_title("%.2e m, %.2e C, %.2e Pa" % (depth, T, P))
    ### second part: plot viscosity vs temperature (e.g. Agard 2016, fig 4)
    # options for temperature
    T_min = 400.0 # C
    T_max = 1000 # C
    T_N = 1000
    # options for viscosity (only affects plotting)
    eta_min = 1e17
    eta_max = 1e25
    # options fot the strain rate
    strain_rates = np.array([1e-15, 1e-13])
    # Options for the pressure
    P = 1.6e3 # 50 km
    Ts = np.linspace(T_min, T_max, T_N)
    for i in range(strain_rates.size):
        strain_rate = strain_rates[i]
        ax = fig.add_subplot(gs[1, i])
        eta_dict = {}  # use for storing all the stress results
        grain_size_dict = {}  # use for stroing all the grain size results
        eta_array_olivine = np.zeros(Ts.shape)
        eta_array_olivine_wet = np.zeros(Ts.shape)
        eta_array_anorthite = np.zeros(Ts.shape)
        eta_array_gabbro_mylonite = np.zeros(Ts.shape)
        eta_array_gabbro_mylonite_small_p = np.zeros(Ts.shape)
        eta_array_gabbro_coarse = np.zeros(Ts.shape)
        eta_array_basalt = np.zeros(Ts.shape)
        for j in range(Ts.size):
            T = Ts[j] + 273.15  # convert to K
            # option of the dry olivine
            eta_diff =  CreepRheology(dry_olivine_diff, strain_rate, P, T, d_olivine_dry, 1.0)
            eta_disl =  CreepRheology(dry_olivine_disl, strain_rate, P, T, d_olivine_dry, 1.0)
            eta_array_olivine[j] = ComputeComposite(eta_diff, eta_disl)
            # option of the wet olivine
            eta_diff = CreepRheologyInAspectViscoPlastic(wet_olivine_diff_aspect, strain_rate, P*1e6, T)
            eta_disl = CreepRheologyInAspectViscoPlastic(wet_olivine_disl_aspect, strain_rate, P*1e6, T)  # use the UI 
            eta_array_olivine_wet[j] = ComputeComposite(eta_diff, eta_disl)
            # for dry anorthite
            d_anorthite_mylonite_dry = d_anorthite_mylonite_dry_array[0]  # use a constant grain size
            eta_diff =  CreepRheology(dry_anorthite_diff, strain_rate, P, T, d_anorthite_mylonite_dry, 1.0)
            eta_disl =  CreepRheology(dry_anorthite_disl, strain_rate, P, T, d_anorthite_mylonite_dry, 1.0)
            eta_array_anorthite[j] = ComputeComposite(eta_diff, eta_disl)
            # gabbro and gabbro with a smaller p
            d_gabbro_mylonite_dry = d_gabbro_mylonite_dry_array[0]  # use a constant grain size
            eta_diff =  CreepRheology(gabbro_diff, strain_rate, P, T, d_gabbro_mylonite_dry, 1.0)
            eta_disl =  CreepRheology(gabbro_disl, strain_rate, P, T, d_gabbro_mylonite_dry, 1.0)
            eta_array_gabbro_mylonite[j] = ComputeComposite(eta_diff, eta_disl)
            eta_diff =  CreepRheology(gabbro_small_p_diff, strain_rate, P, T, d_gabbro_mylonite_dry, 1.0)
            eta_disl =  CreepRheology(gabbro_small_p_disl, strain_rate, P, T, d_gabbro_mylonite_dry, 1.0)
            eta_array_gabbro_mylonite_small_p[j] = ComputeComposite(eta_diff, eta_disl)
            # coarse gabbro
            eta_diff =  CreepRheology(gabbro_diff, strain_rate, P, T, d_gabbro_coarse, 1.0)
            eta_disl =  CreepRheology(gabbro_disl, strain_rate, P, T, d_gabbro_coarse, 1.0)
            eta_array_gabbro_coarse[j] = ComputeComposite(eta_diff, eta_disl)
            # basalt, this only has dislocation creep
            eta_disl =  CreepRheology(basalt_disl, strain_rate, P, T, d_basalt_coarse, 1.0)
            eta_array_basalt[j] = eta_disl 
        eta_dict['olivine_dry'] = eta_array_olivine
        eta_dict['olivine_wet'] = eta_array_olivine_wet
        eta_dict['anorthite_mylonite_dry'] = eta_array_anorthite
        eta_dict['gabbro_mylonite'] = eta_array_gabbro_mylonite
        eta_dict['gabbro_mylonite_small_p'] = eta_array_gabbro_mylonite_small_p
        eta_dict['gabbro_course'] = eta_array_gabbro_coarse
        eta_dict['basalt'] = eta_array_basalt
        grain_size_dict['olivine_dry'] = d_olivine_dry
        grain_size_dict['olivine_wet'] = d_olivine_wet
        grain_size_dict['anorthite_mylonite_dry'] = d_anorthite_mylonite_dry
        grain_size_dict['gabbro_mylonite'] = d_gabbro_mylonite_dry
        grain_size_dict['gabbro_mylonite_small_p'] = d_gabbro_mylonite_dry
        grain_size_dict['gabbro_course'] = d_gabbro_coarse 
        grain_size_dict['basalt'] = d_basalt_coarse
        for key, value in eta_dict.items():
            _label = key + ", d = " + str(grain_size_dict[key]) + " um"
            ax.semilogy(Ts, value, label=_label)
        ax.legend()
        ax.set_xlabel("Temperature (C)")
        ax.set_xlim([T_min, T_max])
        ax.set_ylim([eta_min, eta_max])
        ax.set_ylabel('Viscosity (Pa s)')
        ax.set_title("%.2e s^-1" % strain_rate)
        ax.grid() 
    if option is not None:
        fig_path = os.path.join(RESULT_DIR, "shear_zone_rheology_summary_%s.png" % option) # path to save the figure 
    else:
        fig_path = os.path.join(RESULT_DIR, "shear_zone_rheology_summary.png") # path to save the figure 
    print("%s: generate figure %s" % (Utilities.func_name(), fig_path))
    fig.savefig(fig_path)


def GetMantleScenario(scenario):
    '''
    A wrapper for the rheology scenario in the mantle
    Inputs:
        scenario - the scenario of mantle rheology to use
    '''
    rheology_aspect = None
    viscosity_profile = None
    Operator = RHEOLOGY_OPR()
    all_scenarios = ("wet_modified_medium", "wet_modified_two_d_subduction")
    if  scenario == "two_d_subduction":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=30e3,\
        dVdiff=-5.5e-6, dVdisl=2.12e-6, dAdiff_ratio=0.33333247873, dAdisl_ratio=1.040297619,\
        jump_lower_mantle=15.0)
    elif  scenario == "wet_modified_medium":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03_wet_mod", jump_lower_mantle=15.0)
    elif  scenario == "wet_modified_modified":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=20e3,\
                                                                    dVdiff=-5.5e-6, dVdisl=-1.2e-6, save_profile=1, debug=True)
    elif  scenario == "HK03_wet":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03", jump_lower_mantle=30.0, use_effective_strain_rate=False)
    elif  scenario == "HK03_wet_F":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03", jump_lower_mantle=30.0)
    elif  scenario == "HK03_wet_with_synthetic_adiabat":
        depths = np.linspace(0.0, 2890e3, 3000)
        Operator.depths = depths.copy()
        rho_ref = 3416.0
        cp = 1250.0
        alpha = 3.1e-5
        g = 10.0
        Ts = 1573.0
        MantleAdiabat = MANTLE_ADIABAT(cp, alpha, g, Ts, approx="constant variables")
        Operator.temperatures = MantleAdiabat.Temperature(depths)
        Operator.pressures = rho_ref * g * depths
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03", jump_lower_mantle=30.0, use_effective_strain_rate=False)
    elif  scenario == "HK03_wet_with_synthetic_adiabat_constant_gradient":
        depths = np.linspace(0.0, 2890e3, 3000)
        Operator.depths = depths.copy()
        rho_ref = 3416.0
        cp = 1250.0
        alpha = 3.1e-5
        g = 10.0
        Ts = 1573.0
        MantleAdiabat = MANTLE_ADIABAT(cp, alpha, g, Ts, approx="constant variables and gradient")
        Operator.temperatures = MantleAdiabat.Temperature(depths)
        Operator.pressures = rho_ref * g * depths
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03", jump_lower_mantle=30.0, use_effective_strain_rate=False)
    elif  scenario == "HK03_dry":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="HK03_dry", jump_lower_mantle=30.0, use_effective_strain_rate=False)
    elif  scenario == "AB17":
        source_dir = os.path.join(ASPECT_LAB_DIR, 'tests', 'integration', 'fixtures', 'test_rheology')
        da_file = os.path.join(source_dir, "depth_average.txt")
        assert(os.path.isfile(da_file))
        # read profile
        Operator.ReadProfile(da_file)
        rheology_aspect, viscosity_profile = Operator.MantleRheology(rheology="AB17", jump_lower_mantle=30.0)
    else:
        raise NotImplementedError('scenario needs to be one of ' + str(all_scenarios))
    return rheology_aspect, viscosity_profile
    


def CompareMantleRheology():
    '''
    compare different scenario of mantle rheology
    '''
    # initiate the canvas
    fig = plt.figure(tight_layout=True, figsize=(12, 10))
    gs = gridspec.GridSpec(2, 2)    

    scenarios = ["HK03_wet", "wet_modified_medium", "wet_modified_modified", "two_d_subduction"]
    # scenarios = ["HK03_wet", "HK03_wet_F", "HK03_dry", "two_d_subduction", "AB17"]
    # scenarios = ["HK03_wet_with_synthetic_adiabat", "HK03_wet_with_synthetic_adiabat_constant_gradient"]
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

    ax1 = fig.add_subplot(gs[0, 0])  # plot of temperature, pressure
    ax1_p = ax1.twiny()
    ax2 = fig.add_subplot(gs[0, 1])  # plot of viscosity
    ax3 = fig.add_subplot(gs[1, 0])  # plot of viscosity in the upper mantle
    ax4 = fig.add_subplot(gs[1, 1])  # plot of viscosity in the upper mantle, 1e-13 strain rate
    i = 0
    for scenario in scenarios:
        _color = colors[i]
        _, viscosity_profile = GetMantleScenario(scenario)
        # get all the variables
        eta = viscosity_profile['composite']
        eta_13 = viscosity_profile['composite_13']
        eta_diff = viscosity_profile['diffusion']
        eta_disl = viscosity_profile['dislocation']
        eta_disl_13 = viscosity_profile['dislocation_13']
        depths = viscosity_profile['depth']
        mask_upper = (depths < 660e3)
        Ps = viscosity_profile['P']
        Ts = viscosity_profile['T']
        # plot temperature and pressure
        label_T = None
        label_P = None
        if i==0:
            ln1_1 = ax1.plot(Ts, depths/1e3, color=_color, label='temperature')
            ln1_2 = ax1_p.plot(Ps/1e9, depths/1e3, '--', color=_color, label='pressure')
            lns1 = ln1_1 + ln1_2
        else:
            ax1.plot(Ts, depths/1e3, color=_color)
            ax1_p.plot(Ps/1e9, depths/1e3, '--', color=_color)
        # plot of viscosity
        ax2.semilogx(eta, depths/1e3, color=_color, label=scenario)
        # plot of viscosity in the upper mantle
        label_diff = None
        label_disl = None
        if i== 0:
            label_diff = "diffusion"
            label_disl = "dislocation"
        ax3.semilogx(eta[mask_upper], depths[mask_upper]/1e3, color=_color)
        ax3.semilogx(eta_diff[mask_upper], depths[mask_upper]/1e3, '--', color=_color, label=label_diff)
        ax3.semilogx(eta_disl[mask_upper], depths[mask_upper]/1e3, '-.', color=_color, label=label_disl)
        # plot of viscosity in the upper mantle, strain rate = 1e-13
        label_disl = None
        if i == 0:
            label_disl = "dislocation"
        ax4.semilogx(eta_13[mask_upper], depths[mask_upper]/1e3, color=_color)
        ax4.semilogx(eta_disl_13[mask_upper], depths[mask_upper]/1e3, '-.', color=_color, label=label_disl)
        i += 1

    # options for plot
    # plot of temperature
    ax1.grid()
    ax1.set_xlim([0, 3500.0])
    ax1.set_ylim([3000.0, 0.0])
    ax1_p.set_xlabel('P (GPa)')
    ax1.set_xlabel('T (K)')
    ax1.set_ylabel('Depth (km)')
    ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes,
            size=20, weight='bold')
    labs = [i.get_label() for i in lns1]
    ax1.legend(lns1, labs)
    # plot of viscosity 
    ax2.grid()
    ax2.set_xlim([1e18, 1e24])
    ax2.set_ylim([3000.0, 0.0])
    ax2.set_xlabel('Viscosity (Pa*s)')
    ax2.set_ylabel('Depth (km)')
    ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes,
            size=20, weight='bold')
    ax2.legend()
    # plot of viscosity in the upper mantle 
    ax3.grid()
    ax3.set_xlim([1e18, 1e24])
    ax3.set_ylim([660.0, 0.0])
    ax3.set_xlabel('Viscosity (Pa*s)')
    ax3.set_ylabel('Depth (km)')
    ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes,
            size=20, weight='bold')
    ax3.legend()
    # plot of viscosity in the upper mantle, strain rate = 1e-13
    ax4.grid()
    ax4.set_xlim([1e18, 1e24])
    ax4.set_ylim([660.0, 0.0])
    ax4.set_xlabel('Viscosity (Pa*s)')
    ax4.set_ylabel('Depth (km)')
    ax4.text(-0.1, 1.1, 'D', transform=ax4.transAxes,
            size=20, weight='bold')
    ax4.legend()
    

    # save figure & json files
    o_dir = os.path.join(ASPECT_LAB_DIR, 'results', 'compare_mantle_rheology')
    if os.path.isdir(o_dir):
        rmtree(o_dir)
    os.mkdir(o_dir)
    figure_path = os.path.join(o_dir, 'compare.png')
    fig.savefig(figure_path)
    print("%s: figure saved %s" % (Utilities.func_name(), figure_path))


def AssignAspectViscoPlasticPhaseRheology(visco_plastic_dict, key, idx, diffusion_creep, dislocation_creep):
    '''
    Inputs:
        visco_plastic_dict(dict): options for the viscoplastic module in the aspect material model
        key: name for the composition
        idx: index for the phase
        diffusion_creep: diffusion creep rheology
        dislocation_creep: dislocation creep rheology
    Return:
        visco_plastic_dict(dict): options for the viscoplastic module in the aspect material model after the change
    '''
    assert((diffusion_creep is not None) or (dislocation_creep is not None))
    # diffusion creep
    if diffusion_creep is not None:
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep, use_effective_strain_rate=True)
        visco_plastic_dict["Prefactors for diffusion creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Prefactors for diffusion creep"], key, idx, diffusion_creep_aspect['A'])
        visco_plastic_dict["Grain size exponents for diffusion creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Grain size exponents for diffusion creep"], key, idx, diffusion_creep_aspect['p'])
        visco_plastic_dict["Activation energies for diffusion creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Activation energies for diffusion creep"], key, idx, diffusion_creep_aspect['E'])
        visco_plastic_dict["Activation volumes for diffusion creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Activation volumes for diffusion creep"], key, idx, diffusion_creep_aspect['V'])
    else:
        visco_plastic_dict["Prefactors for diffusion creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Prefactors for diffusion creep"], key, idx, 1e-31)
    # dislocation creep
    if dislocation_creep is not None:
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=True)
        visco_plastic_dict["Prefactors for dislocation creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Prefactors for dislocation creep"], key, idx, dislocation_creep_aspect['A'])
        visco_plastic_dict["Stress exponents for dislocation creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Stress exponents for dislocation creep"], key, idx, dislocation_creep_aspect['n'])
        visco_plastic_dict["Activation energies for dislocation creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Activation energies for dislocation creep"], key, idx, dislocation_creep_aspect['E'])
        visco_plastic_dict["Activation volumes for dislocation creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Activation volumes for dislocation creep"], key, idx, dislocation_creep_aspect['V'])
    else:
        visco_plastic_dict["Prefactors for dislocation creep"] = \
            ReplacePhaseOption(visco_plastic_dict["Prefactors for dislocation creep"], key, idx, 1e-31)
    return visco_plastic_dict



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
    parser.add_argument('-th', '--thermal', type=str,
                        default="hpc",
                        help='the thermal profile to use')
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
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep, use_effective_strain_rate=arg.use_effective_strain_rate)
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=arg.use_effective_strain_rate)
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
        if diffusion_creep is not None:
            eta_diff = CreepRheology(diffusion_creep, arg.strain_rate, arg.pressure, arg.temperature, use_effective_strain_rate=arg.use_effective_strain_rate)
            stress_diff = CreepStress(diffusion_creep, arg.strain_rate, arg.pressure, arg.temperature, 1e4, 1000.0)
            print("eta_diff = %4e" % eta_diff) # screen output
            print("stress_diff = %.4e" % stress_diff)
        if dislocation_creep is not None:
            eta_disl = CreepRheology(dislocation_creep, arg.strain_rate, arg.pressure, arg.temperature, use_effective_strain_rate=arg.use_effective_strain_rate)
            stress_disl = CreepStress(dislocation_creep, arg.strain_rate, arg.pressure, arg.temperature, 1e4, 1000.0)
            print("eta_disl = %4e" % eta_disl) # screen output
            print("stress_disl = %.4e" % stress_disl)
        if (diffusion_creep is not None) and (dislocation_creep is not None):
            eta_diff = CreepRheology(diffusion_creep, arg.strain_rate, arg.pressure, arg.temperature)
            eta_disl = CreepRheology(dislocation_creep, arg.strain_rate, arg.pressure, arg.temperature)
            eta_comp = ComputeComposite(eta_disl, eta_diff)
            print("eta_comp = %4e" % eta_comp) # screen output

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
        DeriveMantleRheology(arg.inputs, save_profile=arg.save_profile, version=arg.version,\
        rheology=arg.rheology, diff=diff, use_effective_strain_rate=arg.use_effective_strain_rate)

    elif _commend == 'plot_strength_profile':
        Operator = STRENGTH_PROFILE(T_type=arg.thermal, brittle="Byerlee", disl="ARCAY17", peierls="MK10")
        # Operator.Execute(creep_type='disl')
        fig_path = os.path.join(ASPECT_LAB_DIR, "results", "strength_profile.png")
        if arg.strain_rate > 0.0:
            PlotStrengh(Operator, fig_path, compute_second_invariant=True, strain_rate=arg.strain_rate)
        else:
            # with an unrealist value of strain_rate, take a range o strain_rate
            PlotStrengh(Operator, fig_path, compute_second_invariant=True)
    
    elif _commend == "plot_shear_zone_rheology_summary":
        ## todo_splot
        PlotShearZoneRheologySummary()
    
    elif _commend == "compare_mantle_rheology":
        CompareMantleRheology()
    elif _commend == "plot_rheology_summary":
        PlotRheologySummaryJson(arg.json)
    else:
        raise CheckValueError('%s is not a valid commend' % _commend)


# run script
if __name__ == '__main__':
    main()