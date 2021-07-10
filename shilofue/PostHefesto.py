# -*- coding: utf-8 -*-
r""" work with Hefesto Phase transtions

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - process hefesto table:

        python -m shilofue.PostHefesto process -i large_data_files/fort.69 -o output/test_table

  - check hefesto table format

        python -m shilofue.PostHefesto check -i large_data_files/fort.69

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


class HEFESTO():

    def __init__(self):
        '''
        initiate class
        '''
        self.header = {}
        self.data = []
        self.version = "1.0.0"
        self.min1 = 0.0 
        self.delta1 = 0.0 
        self.number1 = 0
        self.min2 = 0.0
        self.delta2 = 0.0
        self.number2 = 0

    def read_table(self, path):
        '''
        read data
        '''
        Plotter = LINEARPLOT('hefesto', {})
        print("Reading Header: %s" % path)
        Plotter.ReadHeader(path)
        print("Reading Data: %s" % path)
        Plotter.ReadData(path)
        self.header = Plotter.header
        self.data = Plotter.data

    def CheckHefesto(self, **kwargs):
        '''
        Checkt the Hefesto lookup table
    
        Inputs:
            kwargs: options
                version: version of this file
        Outputs:
            Output to sceen whether the contents are a rectangular table
        '''
        # read dimension info
        print("Read information of the 1st dimension")
        col_P = self.header['Pi']['col']
        min1, delta1, number1 = ReadFirstDimension(self.data[:, col_P])
        print("Dimention 1 has %d entries" % number1)
        print("Checking data")
        is_correct = CheckDataDimension(self.data[:, col_P], min1, delta1, number1)
        if is_correct:
            print('Everything is all right of this file')
        else:
            raise Exception('Something is wrong')
    
    def ProcessHefesto(self, o_path, **kwargs):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            o_path: a output path
            kwargs: options
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        # read dimension info
        col_P = self.header['Pi']['col']
        col_T = self.header['Ti']['col']
        self.min1, self.delta1, self.number1 = ReadFirstDimension(self.data[:, col_P])
        self.min2, self.delta2, self.number2 = ReadSecondDimension(self.data[:, col_T])
        # output
        self.OutputHefesto(o_path)

    def OutputHefesto(self, o_path):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            odata: data to be outputed
            o_path: a output path
            kwargs: options
                version: version of this file
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        print("Outputing Data: %s" % o_path)
        col_P = self.header['Pi']['col']
        col_T = self.header['Ti']['col']
        col_alpagg = self.header['alpagg']['col']
        columns = [col_P, col_T, col_alpagg]
        with open(o_path, 'a') as fout: 
            fout.write(self.version + '\n')  # version
            fout.write(os.path.basename(o_path) + '\n') # filenamea
            fout.write('2\n')  # dimension
            fout.write('P(bar)\n')
            fout.write('\t%s\n' % self.min1)
            fout.write('\t%s\n' % self.delta1)
            fout.write('\t%s\n' % self.number1)
            fout.write('T(K)\n')
            fout.write('\t%s\n' % self.min2)
            fout.write('\t%s\n' % self.delta2)
            fout.write('\t%s\n' % self.number2)
            fout.write('\t%s\n' % len(columns))
            fout.write('%-20s%-20s%-20s\n' % ('P(bar)', 'T(k)', 'alpha,1/K'))  # column info
            np.savetxt(fout, self.data[:, columns], fmt='%-20.8e')
        print("New file generated: %s" % o_path) 

    def PlotHefesto(self):
        '''
        Plot the Hefesto lookup table
    
        Inputs:
            -
        Returns:
            -
        '''
        pass


def CheckDataDimension(nddata, min1, delta1, number1):
    tolerance = 1e-6
    i = 0
    i1 = 0
    is_correct = True
    while True:
        value1 = min1 + delta1 * i1
        if i >= nddata.shape[0]:
            if i1 < number1 - 1:
                print('entry %d(index of row in the data part) is incorrect(missing at the end), the correct value is %.7f' % (i, value1))
            break
        elif i1 >= number1:
            # index in the first dimension exceed maximum
            # move to the next value in the second dimension
            i1 = 0
        elif abs(nddata[i] - value1) > tolerance:
            # value in the first dimension doesn't match
            # shoot message and move over to the next value of the second dimension in our data
            print('entry %d(index of row in the data part) is incorrect(%.7f), the correct value is %.7f' % (i, nddata[i], value1))
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