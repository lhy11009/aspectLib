# -*- coding: utf-8 -*-
r"""Compute the activation volume for lower mantle

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m shilofue.TwoDSubduction0.lower_mantle_V compute

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
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def LowerMantleV():
    '''
    Compute lower mantle activation volume from numbers of lower mantle parameters
    Inputs:
        -
    Returns:
        -
    '''
    E = 3.1700e5
    grad_Tad = (3429.911494552081 - 1971.4774459856458) / (2888.5696990732145 - 746.9661472859572) * 1e-3
    rho = 3500.0
    g = 10.0
    T = 2500.0
    P = 3e10
    V = E * grad_Tad / (rho * g * T - P * grad_Tad)
    print(V)
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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'compute':
        # example:
        LowerMantleV()

# run script
if __name__ == '__main__':
    main()