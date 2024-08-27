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
# from matplotlib import cm
from matplotlib import pyplot as plt

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
  - default usage: \n\
\n\
        python -m \
        ")

class HAGER_CONRAD1999():
    '''
    Analytical solution from Hager and Conrad 1999
    '''
    def __init__(self, eta_l):
        '''
        eta_l - effective slab viscosity
        '''
        self.eta_l = eta_l
        pass

    def ComputeConvergence(self, L_l, h_l, R_l, eta_sz, zeta_f):
        '''
        alpha - thermal expansivity
        rho - density of the slab
        g - gravity acceleration
        Delta_T - temperature difference from slab bottom to top
        L_l - the slab length
        h_l - the slab thickness
        R_l - the slab bending curvature
        zeta - aspect ratio of mantle circulation (cell width / heigth)
        zeta_f - the interface aspect ratio
        eta_m - reference mantle viscosity
        eta_sz - shear zone viscosity
        '''
        rho = 3200.0
        g = 9.8
        alpha = 3e-5
        Delta_T = 800 # K
        B = rho * g * alpha * Delta_T * L_l * h_l
        eta_m = 2.5e20 # Pa s
        zeta = 2.0
        # coefficents defined in Conrad & Hager 1999
        C_s = 1 / np.pi**0.5
        C_l = 2.5
        C_m = 2.5
        C_f = 1.2
        C_zeta = 3 * (zeta + C_m)
        # nondimentionalize the viscosities
        eta_l_dot = self.eta_l / eta_m
        eta_sz_dot = eta_sz / eta_m
        # nondimentionalize the bending radius
        r = R_l / h_l
        Vc = B / eta_m * C_s / (C_zeta + C_l * eta_l_dot * r**(-3.0) + C_f * eta_sz_dot * zeta_f)
        return Vc


class CAP13():
    '''
    class for the Capitanio 13 paper
    '''
    def __init__(self, **kwargs):
        # default values
        self.rho_0 = kwargs.get('rho_0', 3000.0) # Reference density
        self.rho_w = kwargs.get('rho_w', 1000.0) # Water density
        self.alpha = kwargs.get('alpha', 1e-5) # Thermal expansivity
        self.kappa = kwargs.get('kappa', 1e-6) # Thermal diffusivity
        self.delta_T = kwargs.get('delta_T', 1300.0) # Temperature difference
        self.g = kwargs.get('g', 9.81) # Gravity Acceleration
        self.eta_M = kwargs.get('eta_M', 10**21.0) # Upper mantle viscosity
        self.eta_L = kwargs.get('eta_L', 1e4*self.eta_M) # Lithospheric viscosity, by ratio
        self.eta_A = kwargs.get('eta_A', 3e-2*self.eta_M) # Asthenosphere viscosity
        self.delta_rho = kwargs.get('delta_rho', 80.0) # density constract between slab and mantle, not given in the paper
        self.h_A = 500e3 # Asthenosphere thickness

    def CalculateCoefficents(self):
        '''
        compute the coefficents for forces used in this solution
        '''
        # buoyancy force
        # plate thickness
        # h = 2.32 * (kappa * t)**0.5
        # length of the slab to slab dept
        # L = h_D * np.arcsin(phi)
        # buoyancy force to slab dimensions
        # F_B = delta_rho * self.g * L * h
        # F_B = c1 * h_D * np.arcsin(phi) * t**0.5
        c1 = 2.32 * self.delta_rho * self.g * self.kappa**0.5
        
        # ridge push force, from analytical solution
        # F_RP = c2 * t
        c2 = self.g * self.delta_rho * self.kappa * (1.0 + 2.0 * self.delta_rho / np.pi / (self.rho_0 - self.rho_w))
        
        # resisting force
        # the drap coefficient, depend on the two lengths of a ellipsoid of W and L
        # K_M = 24 * 2**0.5 * (1 + np.log(W/L))
        # the resisting force, this includes the shear drag and the anchoring force as the two components
        # F_MD = -c3 * (1 + np.log(W/L)) * U
        c3 = 24 * 2**0.5 * self.eta_M
        
        # the athenospher contour flow resistance
        # tau_A = -2 * eta_A * u / h_A * (2.0 + 3.0 * h / h_A)
        # F_AD = tau_A * l
        # F_AD = -c4 * u * l
        # u: velocity of the plate; l: length of the plate
        # c4 = 2 * eta_A / h_A * (2.0 + 3.0 * h / h_A)
        # assume average values for h / h_A = 0.2
        c4 = 2 * self.eta_A / self.h_A * (2.0 + 3.0 / 5.0)
        
        # resistance to bending
        # minimum bending radius is
        # R = 0.5 * h**3.0 (eta_L / eta_M)**0.5
        # the bending force to radius is
        # F_BR = - 2.0 / 3.0 * U * eta_L * (h / R) **3.0
        # the beding force to velocity is
        # F_BR = - 16.0 / 3.0 * (eta_M / eta_L)**0.5 * eta_M * U
        # assume constant ratio between lithosphere and mantle
        # note here the assumption is eta_M / eta_L = 1.0; Question: is this right?
        # F_BR = -c5 * U
        c5 = 16.0 * self.eta_M / 3.0

        # assign this to a dict
        coeffs = {'c1': c1, 'c2': c2, 'c3': c3, 'c4': c4, 'c5': c5}

        return coeffs
    
    def CalculateVelocities(self, l, t, phi, W, h_D):
        '''
        Calculate the sinking rate, plate sliding, trench retreat velocity and the critical age
        Inputs:
            (all in U/I)
            phi: dipping angle
            l: plate length
            t: plate age
            W: slab width
            L: slab length
            h_D: slab depth
        '''
        # first get the model coefficents
        coeffs = self.CalculateCoefficents()
        c1 = coeffs['c1']
        c2 = coeffs['c2']
        c3 = coeffs['c3']
        c4 = coeffs['c4']
        c5 = coeffs['c5']

        # slab buoyancy is balanced by viscous resistance and bending resistance
        # F_B + F_MD + F_BR = 0.0 -> sinking rate
        L = h_D / np.sin(phi) # length of the slab
        v = c1 * h_D * np.arcsin(phi) * t**0.5 / (c5 - c3 * (1 + np.log(W/L)))

        # ridge push is balanced by the athenospheric drag -> sliding velocity
        # l: length of the plate
        u = c2 * t / c4 / l

        # from the geometric relation (in their eq 2 and 5)
        u_T_A = c1 * h_D * np.cos(phi) * np.sin(phi)**(-2.0) * t**0.5
        u_T_B = c5 - c3 * (1 + np.log(W/L))
        u_T_C = c2 * t / c4 / l
        u_T =  u_T_A / u_T_B - u_T_C

        # a critical age, between trench advance and retreat (eq 6)
        t_c_A = c4 * l / c2
        t_c_B = (c1 * h_D * np.cos(phi) * np.sin(phi)**(-2.0))
        t_c_C = (c5 - c3 * (1.0 + np.log(W/L)))
        t_c = (t_c_A * t_c_B / t_c_C)**2.0

        results = {'v': v, 'u': u, 'u_T': u_T, 't_c': t_c}

        return results, coeffs


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
    elif _commend == 'foo':
        # example:
        SomeFunction('foo')
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()