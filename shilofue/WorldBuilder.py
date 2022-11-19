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

def slab_surface_profile(p0_in, slab_lengths_in, slab_dips_in, coordinate_system, **kwargs):
    '''
    descriptions
    Inputs:
        (note: all the angles are passed in with degree, following the World Builder)
        p0 - a start point, in cartesion or spherical coordinates (radian)
        slab_lengths - segment lengths
        slab_dips - segment dip angles
        coordinate_system: "cartesian" or "spherical"
    Returns:
        p1 - an end point of the segment, in cartesion or spherical coordinates (radian)
    '''
    assert(p0_in.size == 2)
    assert(len(slab_lengths_in) == len(slab_dips_in))
    assert(coordinate_system in ["cartesian", "spherical"])
    # num = len(slab_lengths) + 1
    num = kwargs.get("num", 20)
    ps = np.zeros((num,2))
    if coordinate_system == "cartesian":
        p0 = p0_in
    elif coordinate_system == "spherical":
        p0 = p0_in
        p0[0] *= np.pi / 180.0
    ps[0, :] = p0
    slab_total_length = 0.0
    slab_accumulate_lengths = [0.0]
    slab_dips_at_input_points = [slab_dips_in[0][0] * np.pi / 180.0]
    for i in range(len(slab_lengths_in)):
        length = slab_lengths_in[i]
        slab_total_length += length
        slab_accumulate_lengths.append(slab_total_length)
        slab_dips_at_input_points.append(slab_dips_in[i][1] * np.pi / 180.0)
    intv_slab_length = slab_total_length / (num - 1)
    slab_accumulate_lengths_interpolated = np.linspace(0.0, slab_total_length, num)
    slab_dips_interpolated = np.interp(slab_accumulate_lengths_interpolated, slab_accumulate_lengths, slab_dips_at_input_points)
    for i in range(1, num):
        slab_dip = slab_dips_interpolated[i-1]
        if coordinate_system == "cartesian":
            x0 = ps[i-1, 0]
            y0 = ps[i-1, 1]
            ps[i, 0] = x0 + intv_slab_length * np.cos(slab_dip)
            ps[i, 1] = y0 - intv_slab_length * np.sin(slab_dip)
        elif coordinate_system == "spherical":
            theta0 = ps[i-1, 0]
            r0 = ps[i-1, 1]
            ps[i, 0] = theta0 + np.arcsin(np.cos(slab_dip) / (1 + r0**2.0/intv_slab_length**2.0 - 2 * r0/intv_slab_length * np.sin(slab_dip))**0.5)
            ps[i, 1] = (intv_slab_length**2.0 + r0**2.0 - 2 * r0 * intv_slab_length * np.sin(slab_dip))**0.5
    return ps


def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
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