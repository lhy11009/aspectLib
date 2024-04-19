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

    def __init__(self, Vt, Vm, Vr, L, lbd, mu, n):
        '''
        initiation
        Attributes:
            Vt - 
            Vm - 
            Vr -
            L -
            lbd
        '''
        self.Vt = Vt
        self.Vm = Vm
        self.Vr = Vr
        self.a = L / 2.0
        self.lbd = lbd
        self.mu = mu
        self.n = n
        # self.As = np.zeros(self.n)
        self.As = 0.001 * np.ones(self.n) * mu / lbd**2.0 * (Vt + Vm) # for test
        pass

    def AssignA(self, As):
        '''
        Assign an initial guess of As
        '''
        self.As = As

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
            yi_n = -i / (self.n-1.0) * self.a
            sum += (self.Fxiyi(x, y, xi, yi) if (i == 0) 
                    else self.Fxiyi(x, y, xi, yi_n) + self.Fxiyi(x, y, xi, yi))\
                  * self.lbd**3.0 * Ai / (12.0 * self.mu)
        return sum
    
    def Vxdz_array(self, x):
        '''
        derive an array of Vxdz over different value of yi
        '''
        Vxdzs = np.zeros(self.n)
        for i in range(self.n):
            y = i * self.a / (self.n - 1.0)
            Vxdzs[i] = self.Vxdz(x, y)
        return Vxdzs
    
    def partial_Vxdz_j_partial_Ai(self, x, i, j):
        '''
        partial derivative of Vxdz at the jth point to Ai
        '''
        temp = (j - i) * (j - i) * self.a * self.a / (self.n-1) / (self.n-1)
        return (x * x - temp) / (x * x + temp) * self.lbd * self.lbd / (12 * self.mu * self.Vr)
    
    def Jocobi(self, x):
        '''
        Jocobian of partial Vxdzs_j partial Ai
        '''
        J = np.zeros([self.n, self.n])
        for i in range(self.n):
            for j in range(self.n):
                J[i, j] = (self.partial_Vxdz_j_partial_Ai(x, i, j) if (j == 0)\
                           else self.partial_Vxdz_j_partial_Ai(x, i, j) + self.partial_Vxdz_j_partial_Ai(x, -i, j))
        return J
    
    def PlotVxdz(self, x, ym, **kwargs):
        '''
        Plot the solution
        '''
        ax = kwargs.get("ax", None)
        N = 1000
        ys = np.linspace(0, ym, N)
        Vxdzs = np.zeros(N)
        for i in range(N):
            y = ys[i]
            Vxdzs[i] = self.Vxdz(x, y)
            Vx_dzs_nd = Vxdzs / self.lbd / self.Vr
        if ax is not None:
            ax.plot(ys / self.a, Vxdzs / self.lbd / self.Vr)


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