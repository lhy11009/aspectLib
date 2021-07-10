# -*- coding: utf-8 -*-
r""" work with Hefesto Phase transtions

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - process hefesto table:

        python -m shilofue.PostHefesto process -i large_data_files/fort.69 -o output/test_table

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
from shilofue.Plot import LINEARPLOT
from scipy import interpolate
from shilofue.Utilities import my_assert

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


def PlotHefesto(path):
    '''
    Plot the Hefesto lookup table

    Inputs:
        -
    Returns:
        -
    '''
    # read data
    Plotter = LINEARPLOT('hefesto', {})
    Plotter.ReadHeader(path)
    Plotter.ReadData(path)
    # get data and interpolate
    print(Plotter.header)
    col_pi = Plotter.header['depth_pi']['col']
    pis = Plotter.data[:, col_pi]
    col_T = Plotter.header['Ti']['col']
    Ts = Plotter.data[:, col_T]
    col_PI = Plotter.header['Pi']['col']
    PIs = Plotter.data[:, col_PI]
    # interpolater = interpolate.RectBivariateSpline(pis, Ts, PIs) 
    pass


def ProcessHefesto(path, o_path, **kwargs):
    '''
    Process the Hefesto lookup table for aspect

    Inputs:
        path: a lookup table from Hefesto
        o_path: a output path
        kwargs: options
            version: version of this file
    Outputs:
        Output of this function is the Perplex file form that could be recognized by aspect
    Returns:
        -
    '''
    # read data
    Plotter = LINEARPLOT('hefesto', {})
    print("Reading Header: %s" % path)
    Plotter.ReadHeader(path)
    print("Reading Data: %s" % path)
    Plotter.ReadData(path)
    col_pi = Plotter.header['depth_pi']['col']
    col_T = Plotter.header['Ti']['col']
    col_P = Plotter.header['Pi']['col']
    col_alpagg = Plotter.header['alpagg']['col']
    # read dimension info
    min1, delta1, number1 = ReadFirstDimension(Plotter.data[:, col_P])
    min2, delta2, number2 = ReadSecondDimension(Plotter.data[:, col_T])
    # output
    version = kwargs.get('version', '1.0.0')
    odata = np.zeros((Plotter.data.shape[0], 3))
    odata[:, 0] = Plotter.data[:, col_P]
    odata[:, 1] = Plotter.data[:, col_T]
    odata[:, 2] = Plotter.data[:, col_alpagg]
    print("Outputing Data: %s" % path)
    with open(o_path, 'a') as fout: 
        fout.write(version + '\n')  # version
        fout.write(os.path.basename(o_path) + '\n') # filenamea
        fout.write('2\n')  # dimension
        fout.write('P(bar)\n')
        fout.write('\t%s\n' % min1)
        fout.write('\t%s\n' % delta1)
        fout.write('\t%s\n' % number1)
        fout.write('T(K)\n')
        fout.write('\t%s\n' % min2)
        fout.write('\t%s\n' % delta2)
        fout.write('\t%s\n' % number2)
        fout.write('\t%s\n' % odata.shape[1])
        fout.write('%-20s%-20s%-20s\n' % ('P(bar)', 'T(k)', 'alpha,1/K'))  # column info
        np.savetxt(fout, odata, fmt='%-20.8e')
    print("New file generated: %s" % o_path) 


def CheckHefesto(path, **kwargs):
    '''
    Checkt the Hefesto lookup table

    Inputs:
        path: a lookup table from Hefesto
        kwargs: options
            version: version of this file
    Outputs:
        Output to sceen whether the contents are a rectangular table
    '''
    # read data
    Plotter = LINEARPLOT('hefesto', {})
    print("Reading Header: %s" % path)
    Plotter.ReadHeader(path)
    print("Reading Data: %s" % path)
    Plotter.ReadData(path)
    col_P = Plotter.header['Pi']['col']
    # read dimension info
    print("Read information of the 1st dimension")
    min1, delta1, number1 = ReadFirstDimension(Plotter.data[:, col_P])
    print("Dimention 1 has %d entries" % number1)
    print("Checking data")
    is_correct = CheckDataDimension(Plotter.data[:, col_P], min1, delta1, number1)
    if is_correct:
        print('Everything is all right of this file')
    else:
        raise Exception('Something is wrong')

def CheckDataDimension(nddata, min1, delta1, number1):
    tolerance = 1e-6
    i = 0
    i1 = 0
    is_correct = True
    while True:
        value1 = min1 + delta1 * i1
        if i >= nddata.shape[0]:
            break
        elif i1 >= number1:
            # index in the first dimension exceed maximum
            # move to the next value in the second dimension
            i1 = 0
        elif abs(nddata[i] - value1) > tolerance:
            # value in the first dimension doesn't match
            # shoot message and move over to the next value of the second dimension in our data
            print('entry %d(index of raw in the data part) is incorrect(%.7f), the correct value is %.7f' % (i, nddata[i], value1))
            is_correct = False
            i1 = 0
            while nddata[i] < nddata[i+1]:
                i += 1
            i += 1
        else:
            # print(i, nddata[i], value1)  # debug
            # move to the next value
            i1 += 1
            i += 1
    return is_correct


def ReadFirstDimension(nddata):
    '''
    Process the data in the fisrt dimension(min, delta, number)
    
    Inputs:
        nddata: a ndarray of the first dimension data
    Returns:
        min: min value
        delta: data interval
        number: number in a column
    '''
    # min value
    min = nddata[0]
    # delta
    delta = nddata[1] - nddata[0]
    # number
    number = nddata.size
    for i in range(0, nddata.size-1):
        if nddata[i] > nddata[i+1]:
            number = i + 1
            break
    return min, delta, number


def ReadSecondDimension(nddata):
    '''
    Process the data in the second dimension(min, delta, number)
    
    Inputs:
        nddata: a ndarray of the first dimension data
    Returns:
        min: min value
        delta: data interval
        number: number in a column
    '''
    # min value
    min = nddata[0]
    # delta
    tolerance = 1e-6
    delta = 0.0
    for i in range(0, nddata.size-1):
        if abs(nddata[i] - nddata[i+1]) / abs(nddata[i]) > tolerance:
            delta = nddata[i+1] - nddata[i]
            sub_size = i + 1
            break
    # number
    print('subsize: ', sub_size, ", size: ", nddata.size)
    my_assert(nddata.size % sub_size == 0, ValueError, 'the table is not regular(rectangle)')
    number = nddata.size // sub_size
    return min, delta, number


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
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='Some outputs')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot':
        # Plot the Hefesto lookup table
        PlotHefesto(arg.inputs)
    
    elif _commend == 'process':
        # Process the Hefesto lookup table for aspect
        ProcessHefesto(arg.inputs, arg.outputs)

    elif _commend == 'check':
        # Checkt the Hefesto lookup table
        CheckHefesto(arg.inputs)

# run script
if __name__ == '__main__':
    main()