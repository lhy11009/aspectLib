# -*- coding: utf-8 -*-
r""" work with Hefesto Phase transtions

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - process hefesto table:

        python -m shilofue.PostHefesto process -i large_data_files/fort.56 -i1 4 -i2 4 -o output/test_table

            i1, and i2 are intervals in the first and second dimension, default is 1.

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

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


class LOOKUP_TABLE():

    def __init__(self):
        '''
        initiate class
        '''
        self.header = {}
        self.data = []
        self.version = "1.0.0"
        self.UnitConvert = Utilities.UNITCONVERT()
        self.min1 = 0.0 
        self.delta1 = 0.0 
        self.number1 = 0
        self.min2 = 0.0
        self.delta2 = 0.0
        self.number2 = 0
        self.indexes = []  # indexes of output data
        self.number_out1 = 0 # number of output
        self.number_out2 = 0
        self.delta_out1 = 0.0  # intervals used to outptu
        self.delta_out2 = 0.0 
        self.oheader = { 'Temperature': 'T(K)',  'Pressure': 'P(bar)' ,  'Density': 'rho,kg/m3',\
        'Thermal_expansivity': 'alpha,1/K', 'Isobaric_heat_capacity': 'cp,J/K/kg',\
        'VP': 'vp,km/s', 'VS': 'vs,km/s', 'Enthalpy': 'h,J/kg' }
        # unit to output
        self.ounit = {'Temperature': 'K', 'Pressure': 'bar', 'Thermal_expansivity': '1/K',\
        'Isobaric_heat_capacity': 'J/K/kg', 'Density': 'kg/m3', 'VP':'km/s', 'VS':'km/s', 'Enthalpy': 'J/kg'}

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

    def Check(self, first_dimension_name, **kwargs):
        '''
        Checkt the Hefesto lookup table, this only check that the first dimension
        is aligned and no data is missings.
    
        Inputs:
            kwargs: options
                version: version of this file
        Outputs:
            Output to sceen whether the contents are a rectangular table
        '''
        # read dimension info
        print("Read information of the 1st dimension")
        col_first = self.header[first_dimension_name]['col']
        min1, delta1, number1 = ReadFirstDimension(self.data[:, col_first])
        print("Dimention 1 has %d entries" % number1)
        print("Checking data")
        is_correct = CheckDataDimension(self.data[:, col_first], min1, delta1, number1)
        if is_correct:
            print('Everything is all right of this file')
        else:
            raise Exception('Something is wrong')
    
    def Process(self, field_names, o_path, **kwargs):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            o_path: a output path
            kwargs: options
                interval1 & 2: interval in the first & second dimension
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        # read dimension info
        first_dimension_name = kwargs.get('first_dimension', 'Pressure')
        second_dimension_name = kwargs.get('second_dimension', 'Temperature')
        col_first = self.header[first_dimension_name]['col']
        col_second = self.header[second_dimension_name]['col']
        self.min1, self.delta1, self.number1 = ReadFirstDimension(self.data[:, col_first])
        self.min2, self.delta2, self.number2 = ReadSecondDimension(self.data[:, col_second])
        # output
        interval1 = kwargs.get('interval1', 1)
        interval2 = kwargs.get('interval2', 1)
        self.indexes = self.IndexesByInterval(interval1, interval2)  # work out indexes
        self.number_out1 = int(np.ceil(self.number1 / interval1)) # number of output
        self.number_out2 = int(np.ceil(self.number2 / interval2))
        # output intervals
        self.delta_out1 = self.delta1 * interval1 # output intervals
        self.delta_out2 = self.delta2 * interval2 # output intervals
        self.OutputHefesto(field_names, o_path)

    def OutputHefesto(self, field_names, o_path):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            odata: data to be outputed
            o_path: a output path
            kwargs: options
                version: version of this file
            field_names: field_name to output, the first two are the first and second dimension
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        UnitConvert = Utilities.UNITCONVERT()
        print("Outputing Data: %s" % o_path)
        # columns
        print("Outputing fields: %s" % field_names)  # debug
        print('first dimension: ', self.number_out1, ", second dimension: ", self.number_out2, ", size:", self.number_out1 * self.number_out2)
        Utilities.my_assert(len(field_names) >= 2, ValueError, 'Entry of field_names must have more than 2 components')
        columns = []

        missing_last = self.data.shape[1]
        missing_fix_values = []
        for field_name in field_names:
            # attach size(field_names) if failed
            try:
                columns.append(self.header[field_name]['col'])
            except KeyError:
                # first check that T or P is not missing
                # then append an imaginary column
                if field_name == 'Temperature':
                    raise KeyError('Abort: Temperature field is missing')
                elif field_name == 'Pressure':
                    raise KeyError('Abort: Pressure field is missing')
                else:
                    # assign an append value
                    print('field %s is missing, going to append manually' % field_name)
                    columns.append(missing_last)
                    missing_last += 1
                    # ask for value
                    missing_fix_value = float(input('Input value:'))
                    missing_fix_values.append(missing_fix_value)
        unit_factors = []
        for field_name in field_names:
            # attach 1 if failed
            try:
                unit_factors.append(self.UnitConvert(self.header[field_name]['unit'], self.ounit[field_name]))
            except KeyError:
                unit_factors.append(1.0)
        # check the output values
        # note that self.indexes[self.number_out1] gives the index of the second member in the second dimension
        tolerance = 1e-5
        temp1 = self.data[self.indexes[1], columns[0]] - self.data[self.indexes[0], columns[0]]
        temp2 = self.data[self.indexes[self.number_out1], columns[1]] - self.data[self.indexes[0], columns[1]]  
        Utilities.my_assert( (abs(temp1 - self.delta_out1) / self.delta_out1) < tolerance,
        ValueError, "Output interval(self.delta_out1) doesn't match the interval in data")
        Utilities.my_assert( (abs(temp2 - self.delta_out2) / self.delta_out2) < tolerance,
        ValueError, "Output interval(self.delta_out2) doesn't match the interval in data")
        # mend self.data if needed
        if missing_last > self.data.shape[1]:
            print("Concatenating missing data")
            new_data = np.ones((self.data.shape[0], missing_last - self.data.shape[1])) *  missing_fix_values
            self.data = np.concatenate((self.data, new_data), axis=1)
        # output
        with open(o_path, 'w') as fout: 
            fout.write(self.version + '\n')  # version
            fout.write(os.path.basename(o_path) + '\n') # filenamea
            fout.write('2\n')  # dimension
            fout.write('%s\n' % self.oheader[field_names[0]])
            fout.write('\t%.8f\n' % (float(self.min1) * unit_factors[0])) # min value
            fout.write('\t%.8f\n' % (float(self.delta_out1) * unit_factors[0]))  # difference, use the output value
            fout.write('\t%s\n' % self.number_out1)  # number of output
            fout.write('%s\n' % self.oheader[field_names[1]])
            fout.write('\t%.8f\n' % (float(self.min2) * unit_factors[1]))
            fout.write('\t%.8f\n' % (float(self.delta_out2) * unit_factors[1]))
            fout.write('\t%s\n' % self.number_out2)
            fout.write('\t%s\n' % len(columns))
            temp = ''
            for field_name in field_names:
                temp += '%-20s' % self.oheader[field_name]
            temp += '\n'
            fout.write(temp)
            # data is indexes, so that only part of the table is output
            np.savetxt(fout, self.data[np.ix_(self.indexes, columns)] * unit_factors, fmt='%-19.8e')
        print("New file generated: %s" % o_path) 

    def IndexesByInterval(self, interval1, interval2):
        '''
        Work out indexes by giving interval(default is 1, i.e. consecutive)
        '''
        Utilities.my_assert(type(interval1) == int and type(interval2) == int, TypeError, "interval1(%s) or interval2(%s) is not int" % (interval1, interval2))
        # indexes in 
        indexes_1 = range(0, self.number1, interval1) 
        indexes_2 = range(0, self.number2, interval2)
        # work out the overall indexes
        indexes = []
        for index_2 in indexes_2:
            for index_1 in indexes_1: 
                indexes.append(index_1 + self.number1 * index_2)
        return indexes

    def PlotHefesto(self):
        '''
        Plot the Hefesto lookup table
    
        Inputs:
            -
        Returns:
            -col_alpagg
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
    Utilities.my_assert(nddata.size % sub_size == 0, ValueError, 'the table is not regular(rectangle)')
    number = nddata.size // sub_size
    print('first dimension: ', sub_size, ", second dimension: ", number, ", total size: ", nddata.size)
    return min, delta, number


def ProcessHefesto(filein, fileout, interval1, interval2):
    # input file
    assert(os.path.isfile(filein))
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.read_table(filein)
    # fields to read in
    field_names = ['Pressure', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.Process(field_names, fileout, interval1=interval1, interval2=interval2)
    # assert something 
    assert(os.path.isfile(fileout))


def CheckHefesto(filein, first_dimension_name):
    '''
    check the data format of a HeFESTo lookup table
    '''
    # input file
    assert(os.path.isfile(filein))
    Hefesto = LOOKUP_TABLE()
    Hefesto.read_table(filein)
    Hefesto.Check(first_dimension_name)


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
    parser.add_argument('-i1', '--interval1', type=int,
                        default=1,
                        help='Interval in the first dimension')
    parser.add_argument('-i2', '--interval2', type=int,
                        default=1,
                        help='Interval in the second dimension')
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
        ProcessHefesto(arg.inputs, arg.outputs, arg.interval1, arg.interval2)

    elif _commend == 'check':
        # Check the Hefesto lookup table
        CheckHefesto(arg.inputs, 'Pressure')

# run script
if __name__ == '__main__':
    main()