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