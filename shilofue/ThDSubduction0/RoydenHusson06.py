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

def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
    pass

class HELE_SHAW():

    def __init__(self, Vt, Vm, L, lbd, mu):
        '''
        initiation
        Attributes:
            Vt - 
            Vm - 
            L -
            lbd
        '''
        self.Vt = Vt
        self.Vm = Vm
        self.a = L / 2.0
        self.lbd = lbd
        self.mu = mu
        self.n = 100
        # self.As = np.zeros(self.n)
        self.As = 0.001 * np.ones(self.n) * mu / lbd**2.0 * (Vt + Vm) # for test
        pass

    def Fxiyi(self, x, y, xi, yi):
        '''
        coefficients of the solution for Vx
        '''
        temp = ((x - xi)**2.0 - (y - yi)**2.0) / ((x - xi)**2.0 + (y - yi)**2.0)
        return temp

    def Vxdz(self, x, y):
        '''
        integration of Vx over z
        '''
        sum = self.lbd * (self.Vt + self.Vm) / 2.0
        for i in range(self.n):
            Ai = self.As[i]
            xi = 0.0
            yi = i / (self.n-1.0) * self.a
            sum += (1.0 if (i == 0) else 2.0) * self.Fxiyi(x, y, xi, yi) * self.lbd**3.0 * Ai / (12.0 * self.mu)
        return sum



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