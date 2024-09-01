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
from operator import index
import numpy as np
import sys, os, argparse
import json, re
# import pathlib
# import subprocess
import numpy as np
from shutil import copy2, rmtree, copytree
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec
from shilofue.Plot import LINEARPLOT
import shilofue.ParsePrm as ParsePrm
from scipy import interpolate

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


class HEFESTO_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with CASE
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("HeFESTo Repository", str, ["HeFESTo repository"], ".", nick='hefesto_dir')
        self.add_key("Output directory", str, ["output directory"], ".", nick='o_dir')
        self.add_key("Case name", str, ["case name"], "foo", nick='case_name')
        self.add_key("Number of processor", int, ["nproc"], 1, nick='nproc')
        self.add_key("Path of a control file", str, ["control path"], ".", nick='control_path')
        self.add_key("dimension T: T1", float, ["T", "T1"], 0.0, nick='T1')
        self.add_key("dimension T: T2", float, ["T", "T2"], 0.0, nick='T2')
        self.add_key("dimension T: nT", int, ["T", "nT"], 0, nick='nT')
        self.add_key("dimension P: P1", float, ["P", "P1"], 0.0, nick='P1')
        self.add_key("dimension P: P2", float, ["P", "P2"], 0.0, nick='P2')
        self.add_key("dimension P: nP", int, ["P", "nP"], 0, nick='nP')
        self.add_key("dimension T: variable", str, ["T", "variable"], 'temperature', nick='T_variable')
        self.add_key("path of the slurm file", str, ["slurm path"], 'foo', nick='slurm_file_base')
        self.add_key("path of the parameter directory", str, ["prm directory"], 'foo', nick='prm_dir')
        self.add_key("split by", str, ["split by"], 'P', nick='split_by')


    def check(self):
        '''
        check inputs are valid
        '''
        hefesto_dir = Utilities.var_subs(self.values[0])
        o_dir = Utilities.var_subs(self.values[1])
        Utilities.my_assert(os.path.isdir(o_dir), FileNotFoundError, "No such directory: %s" % o_dir)
        control_path = Utilities.var_subs(self.values[4])
        Utilities.my_assert(os.path.isfile(control_path), FileNotFoundError, "No such file: %s" % control_path)
        # order of T1, T2; P1, P2
        T1 = self.values[5]
        T2 = self.values[6]
        assert(T1 <= T2)
        P1 = self.values[8]
        P2 = self.values[9]
        assert(P1 <= P2)
        T_variable = self.values[11]
        assert(T_variable in ["temperature", "entropy"])
        slurm_file_base = Utilities.var_subs(self.values[12])
        assert(os.path.isfile(slurm_file_base))


    def to_distribute_parallel_control(self):
        '''
        Interface to the DistributeParallelControl function
        '''
        hefesto_dir = Utilities.var_subs(self.values[0])
        o_dir = Utilities.var_subs(self.values[1])
        case_name = Utilities.var_subs(self.values[2])
        nproc = self.values[3]
        control_path = Utilities.var_subs(self.values[4])
        T1 = self.values[5]
        T2 = self.values[6]
        nT = self.values[7]
        P1 = self.values[8]
        P2 = self.values[9]
        nP = self.values[10]
        T_variable = self.values[11]
        slurm_file_base = Utilities.var_subs(self.values[12])
        prm_dir = Utilities.var_subs(self.values[13])
        split_by = Utilities.var_subs(self.values[14])
        return hefesto_dir, o_dir, case_name, nproc, control_path,\
            T1, T2, nT, P1, P2, nP, T_variable, slurm_file_base, prm_dir,\
            split_by

    def case_dir(self):
        '''
        return the path to the case directory
        '''
        o_dir = Utilities.var_subs(self.values[1])
        case_name = Utilities.var_subs(self.values[2])
        case_dir = os.path.join(o_dir, case_name)
        return case_dir


class LOOKUP_TABLE():

    def __init__(self):
        '''
        initiate class
        '''
        self.header = {}
        self.data = []
        self.version = "1.0.0"
        self.UnitConvert = Utilities.UNITCONVERT()
        self.first_dimension_name = None
        self.min1 = 0.0 
        self.delta1 = 0.0 
        self.number1 = 0
        self.second_dimension_name = None
        self.min2 = 0.0
        self.delta2 = 0.0
        self.number2 = 0
        self.indexes = []  # indexes of output data
        self.number_out1 = 0 # number of output
        self.number_out2 = 0
        self.delta_out1 = 0.0  # intervals used to outptu
        self.delta_out2 = 0.0 
        self.fort_56_header = { "Pressure": {"col": 0, "unit": "GPa"}, "Depth": {"col": 1, "unit": "km"},\
                               "Temperature": {"col": 2, "unit": "K"}, "Density": {"col": 3, "unit": "g/cm^3"},\
                                "Bulk sound Velocity": {"col": 4, "unit": "km/s"}, "VS": {"col": 5, "unit": "km/s"},\
                                "VP": {"col": 6, "unit": 'km/s'}, "VS modified by attenuation": {"col": 7, "unit": 'km/s'},\
                                "VP modified by attenuation": {"col": 8, "unit": 'km/s'}, "Enthalpy":{"col": 9, "unit": 'kJ/g'},\
                                "Entropy": {"col": 10, "unit": "J/g/K"}, "Thermal_expansivity": {"col": 11, "unit": "10-5 K-1"},\
                                "Isobaric_heat_capacity": {"col": 12, "unit": "J/g/K"}, "isothermal bulk modulus": {"col": 13, "unit": "GPa"},\
                                "Shear attenuation": {"col": 14, "unit": "1"}, "Longitudinal attenuation":{"col": 15, "unit": "1"},\
                                "Quenched density": {"col": 16, "unit": "g/cm^3"}, "most abundant phase": {"col": 17, "unit": None}}
        self.oheader = { 'Temperature': 'T(K)',  'Pressure': 'P(bar)' ,  'Density': 'rho,kg/m3',\
        'Thermal_expansivity': 'alpha,1/K', 'Isobaric_heat_capacity': 'cp,J/K/kg',\
        'VP': 'vp,km/s', 'VS': 'vs,km/s', 'Enthalpy': 'h,J/kg', "Entropy": 's,J/K/kg',\
        'Omph_GHP': 'Omph', "Gt_HGP": "Gt"}
        # header for perplex table
        self.perplex_header = { 'T': "Temperature", "P": 'Pressure', "rho": 'Density',\
                               "alpha": 'Thermal_expansivity', 'cp': "Isobaric_heat_capacity",\
                                'vp': "VP", 'vs': "VS", 'h': "Enthalpy", "s": "Entropy", "H2O": "H2O",\
                                "cAmph(G)": "cAmph_G", "Ep(HP11)": "Ep_HP11", "law": "law", "Omph(GHP)": "Omph_GHP", "Gt(HGP)": "Gt_HGP"}
        # unit to output
        self.ounit = {'Temperature': 'K', 'Pressure': 'bar', 'Thermal_expansivity': '1/K',\
        'Isobaric_heat_capacity': 'J/K/kg', 'Density': 'kg/m3', 'VP':'km/s', 'VS':'km/s', 'Enthalpy': 'J/kg', "Entropy": "J/K/kg"}
        # pressure entropy table
        self.Pinterp_rows = None
        self.PS_data = None
        self.PS_number_out1 = 0 # number of output
        self.PS_number_out2 = 0
        self.PS_delta_out1 = 0.0  # intervals used to outptu
        self.PS_delta_out2 = 0.0 


    def ReadHeFestoTable(self, path):
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


    def ReadPerplex(self, path, **kwargs):
        '''
        read Perplex data
        Inputs:
            path(str):
                file path of the table
            kwargs:
                header_rows : number of header raws
        '''
        header_rows = kwargs.get('header_rows', 0)
        line = ""
        with open(path, 'r') as fin:
            for i in range(header_rows):
                line = fin.readline()
                assert(line != "")
        header, unit = ParsePerplexHeader(line)

        # match the header
        for i in range(len(header)):
            header_to = self.perplex_header[header[i]]
            self.header[header_to] = {}
            self.header[header_to]['col'] = i
            self.header[header_to]['unit'] = unit[i]
        
        # read data
        self.data = np.loadtxt(path, skiprows=header_rows) 

    def ReadRawFort56(self, _path):
        '''
        Read a raw fort.56 file
        Inputs:
            _path: _path to the input fort.56
        '''
        print("%s: reading file" % Utilities.func_name())
        assert(os.path.isfile(_path))  # assert file exists
        self.header = self.fort_56_header
        self.data = np.genfromtxt(_path)
        print("%s: data dimension: " % Utilities.func_name(), self.data.shape)

    def AllFields(self):
        '''
        return all fields options
        '''
        fields = []
        for key, _ in self.header.items():
            fields.append(key)
        return fields
    
    def export_fort56_vss(self):
        '''
        export the vss data from the fort.56 file
        '''
        col_depth = self.header["Depth"]["col"]
        col_vs = self.header["VS"]["col"]
        depths = self.data[:, col_depth]
        Vss = self.data[:, col_vs]
        return depths, Vss
    
    def export_density_profile(self):
        '''
        export density profile
        '''
        col_depth = self.header["Depth"]["col"]
        col_density = self.header["Density"]["col"]
        depths = self.data[:, col_depth]
        densities = self.data[:, col_density]
        return depths, densities
    
    def export_temperature_profile(self):
        '''
        export temperature profile
        '''
        col_depth = self.header["Depth"]["col"]
        col_T = self.header["Temperature"]["col"]
        depths = self.data[:, col_depth]
        temperatures = self.data[:, col_T]
        return depths, temperatures
    
    def export_thermal_expansivity_profile(self):
        '''
        export density profile
        '''
        col_depth = self.header["Depth"]["col"]
        col_alpha = self.header["Thermal_expansivity"]["col"]
        depths = self.data[:, col_depth]
        alphas = self.data[:, col_alpha]
        return depths, alphas

    def export_pressure_profile(self):
        '''
        export density profile
        '''
        col_depth = self.header["Depth"]["col"]
        col_pressure = self.header["Pressure"]["col"]
        depths = self.data[:, col_depth]
        pressures = self.data[:, col_pressure]
        return depths, pressures
    
    def export_field_profile(self, field):
        '''
        export field profile
        Inputs:
            field (str): the field to output
        '''
        col_depth = self.header["Depth"]["col"]
        col_field = self.header[field]["col"]
        depths = self.data[:, col_depth]
        _data = self.data[:, col_field]
        return depths, _data

    def export_field(self, field):
        '''
        export field
        Inputs:
            field (str): the field to output
        '''
        col_field = self.header[field]["col"]
        _data = self.data[:, col_field]
        return _data
    
    def export_field_mesh(self, field):
        '''
        export field in a 2-d mesh
        '''
        col_field = self.header[field]["col"]

        # call update and process with information of the mesh
        self.Update()
        col_first = self.header[self.first_dimension_name]['col']
        col_second = self.header[self.second_dimension_name]['col']

        # initiate 2-d arrays 
        D1 = np.zeros((self.number1, self.number2))
        D2 = np.zeros((self.number1, self.number2))
        F = np.zeros((self.number1, self.number2))

        # process data
        for i in range(self.number1):
            for j in range(self.number2):
                D1[i, j] = self.data[i, col_first]
                D2[i, j] = self.data[j*self.number1, col_second]
                F[i, j] = self.data[j*self.number1 + i, col_field]
        
        return D1, D2, F
    
    def fix_field_nan_value(self, field, value):
        '''
        fix nan value of field
        Inputs:
            field (str): name of field
            value (float): value to fix with
        '''
        col_field = self.header[field]["col"]
        for i in range(self.data.shape[0]):
            # if self.data[i, col_field] == float('nan'):
            if np.isnan(self.data[i, col_field]):
                self.data[i, col_field] = value


    def Update(self, **kwargs):
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

        # get first dimension name
        find_first = False
        find_second = False
        for key, value in self.header.items():
            if value['col'] == 0:
                self.first_dimension_name = key
                find_first = True
            elif value['col'] == 1:
                self.second_dimension_name = key
                find_second = True
        assert(find_first and find_second)

        # update and check the first dimension
        col_first = self.header[self.first_dimension_name]['col']
        min1, delta1, number1 = ReadFirstDimension(self.data[:, col_first])
        self.min1 = min1
        self.delta1 = delta1
        self.number1 = number1
        print("Dimention 1 has %d entries" % number1)
        print("Checking data")
        is_correct = CheckDataDimension(self.data[:, col_first], min1, delta1, number1)
        if is_correct:
            print('Everything is all right of this file')
        else:
            raise Exception('Something is wrong')
        
        # update the second dimension information
        col_second = self.header[self.second_dimension_name]['col']
        min2, delta2, number2 = ReadSecondDimension(self.data[:, col_second])
        self.min2 = min2
        self.delta2 = delta2
        self.number2 = number2

    def CreateNew(self, new_data, _name, oheader):
        # todo_new
        self.header[_name] = {'col': self.data.shape[1]}
        self.oheader[_name] = oheader
        self.ounit[_name] = oheader
        self.data = np.concatenate((self.data, new_data), axis=1)
    
    def Process(self, field_names, o_path, **kwargs):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            o_path: a output path
            kwargs: options
                interval1 & 2: interval in the first & second dimension
                digit: digit of output numbers
                file_type: type of the output file, perple_x or structured
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        first_dimension_name = kwargs.get('first_dimension', 'Pressure')
        second_dimension_name = kwargs.get('second_dimension', 'Temperature')
        fix_coordinate_minor = kwargs.get("fix_coordinate_minor", False)
        interval1 = kwargs.get('interval1', 1)
        interval2 = kwargs.get('interval2', 1)

        # read dimension info
        col_first = self.header[first_dimension_name]['col']
        col_second = self.header[second_dimension_name]['col']
        self.min1, self.delta1, self.number1 = ReadFirstDimension(self.data[:, col_first])
        self.min2, self.delta2, self.number2 = ReadSecondDimension(self.data[:, col_second])

        # fix minor error in the coordinate data
        if fix_coordinate_minor:
            FixFirstDimensionMinor(self.data[:, col_first], self.number1)
            FixSecondDimensionMinor(self.data[:, col_second], self.number1)

        # output
        if interval1 == 1 and interval2 == 1:
            self.indexes = np.array(range(self.number1*self.number2))
        else:
            print("%s begin indexing" % Utilities.func_name())  # debug
            self.indexes = np.array(self.IndexesByInterval(interval1, interval2))  # work out indexes
            print("%s finish indexing" % Utilities.func_name())  # debug
        self.number_out1 = int(np.ceil(self.number1 / interval1)) # number of output
        self.number_out2 = int(np.ceil(self.number2 / interval2))
        # output intervals
        self.delta_out1 = self.delta1 * interval1 # output intervals
        self.delta_out2 = self.delta2 * interval2 # output intervals
        self.OutputPerplexTable(field_names, o_path, **kwargs)

    def OutputPerplexTable(self, field_names, o_path, **kwargs):
        '''
        Process the Hefesto lookup table for aspect
    
        Inputs:
            o_path: a output path
            field_names: field_name to output, the first two are the first and second dimension
            kwargs: options
                version: version of this file
                digit: digit of output numbers
                file_type: type of the output file, perple_x or structured
                exchange_dimension: exchange the 1st and 2nd dimensions
        Outputs:
            Output of this function is the Perplex file form that could be recognized by aspect
        Returns:
            -
        '''
        digit = kwargs.get('digit', 8)
        file_type = kwargs.get("file_type", 'perple_x')
        exchange_dimension = kwargs.get("exchange_dimension", False)
        assert(file_type in ['perple_x', 'structured'])

        UnitConvert = Utilities.UNITCONVERT()
        print("%s: Outputing Data" % Utilities.func_name())
        # columns
        print("Outputing fields: %s" % field_names)
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
            if file_type == "perple_x":
                if exchange_dimension:
                    raise NotImplementedError()
                # write header for perple_x file
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
            elif file_type == "structured":
                # write header for structured file
                fout.write("# This is a data output from HeFESTo\n")
                if exchange_dimension:
                    fout.write("# Independent variables are %s and %s\n" % (self.oheader[field_names[1]], self.oheader[field_names[0]]))
                    fout.write("# POINTS: %d %d\n" % (self.number_out2, self.number_out1)) 
                else:
                    fout.write("# Independent variables are %s and %s\n" % (self.oheader[field_names[0]], self.oheader[field_names[1]]))
                    fout.write("# POINTS: %d %d\n" % (self.number_out1, self.number_out2)) 
            temp = ''
            if exchange_dimension:
                temp += '%-20s%-20s' % (self.oheader[field_names[1]], self.oheader[field_names[0]])
            else:
                temp += '%-20s%-20s' % (self.oheader[field_names[0]], self.oheader[field_names[1]])
            for i in range(2, len(field_names)):
                field_name = field_names[i]
                temp += '%-20s' % self.oheader[field_name]
            temp += '\n'
            fout.write(temp)
            # data is indexes, so that only part of the table is output
            _format = '%-' + str(digit + 11) + '.' + str(digit) + 'e'
            indexes_output = None; columns_output = None; unit_factors_output = None
            if exchange_dimension:
                indexes_output = ExchangeDimensions(self.indexes, self.number_out1, self.number_out2)
                columns_output = columns.copy()
                columns_output[0] = columns[1]
                columns_output[1] = columns[0]
                unit_factors_output = unit_factors.copy()
                tmp = unit_factors[0]
                unit_factors_output[0] = unit_factors[1]
                unit_factors_output[1] = tmp
            else:
                indexes_output = self.indexes
                columns_output = columns
                unit_factors_output = unit_factors
            np.savetxt(fout, self.data[np.ix_(indexes_output, columns_output)] * unit_factors_output, fmt=_format)
        print("New file generated: %s" % o_path) 
    
    def OutputPressureEntropyTable(self, field_names, o_path):
        '''
        Converts the dataset and generate a lookup table
        from (P, S) -> T
            o_path: a output path
            field_names: field_name to output, the first two are the first and second dimension
        '''
        # initiation
        UnitConvert = Utilities.UNITCONVERT()
        Utilities.my_assert(len(field_names) == self.PS_data.shape[1], ValueError,\
                            'Entry of field_names must equal the size of the pressure-entropy dataset')

        # figure out the factors of unit converting 
        unit_factors = []
        for field_name in field_names:
            # attach 1 if failed
            try:
                unit_factors.append(self.UnitConvert(self.header[field_name]['unit'], self.ounit[field_name]))
            except KeyError:
                unit_factors.append(1.0)

        # write output
        with open(o_path, 'w') as fout:
            fout.write("# This is a data output from.\n")
            fout.write("# Independent variables are entropy and pressure.\n")
            fout.write("# POINTS: %s %s\n" % (self.PS_number_out1, self.PS_number_out2))
            temp = ''
            for field_name in field_names:
                temp += '%-20s' % self.oheader[field_name]
            temp += '\n'
            fout.write(temp)
            np.savetxt(fout, self.PS_data * unit_factors, fmt='%-19.8e')
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
        i2 = 0 # used for printing the percentage of completeness
        last_ratio = 0.0
        for index_2 in indexes_2:
            ratio = i2 / len(indexes_2)
            if ratio > last_ratio + 0.01:
                last_ratio = ratio
                print("percent = %.2f" % (ratio*100), end='\r')
            for index_1 in indexes_1: 
                indexes.append(index_1 + self.number1 * index_2)
            i2 += 1
        return indexes

    def InterpolatePressureEntropyByIndex(self, index_p, entropies, field_names, PS_rows, **kwargs):
        '''
        Interpolate the data with a pressure index
        Inputs:
            index_p (int): pressure index
            entropies (list): entropy inputs
            field_names (list): names of field to interpolate
            PS_rows (list): range of rows to enter into the PS_data 2-d ndarray
        kwargs:
            debug: run debug mode
        '''
        # initiate
        assert(self.PS_data is not None) # PS_data is initiated
        col_entropy = self.header["Entropy"]['col']
        col_pressure = self.header["Pressure"]['col']
        pressure = 0.0
        debug = kwargs.get("debug", False)

        # row index for this pressure
        if self.first_dimension_name == "Pressure":
            if self.Pinterp_rows is None:
                self.Pinterp_rows = np.zeros(self.number2).astype(int)
            for i in range(self.number2):
                row = self.number1 * i + index_p 
                self.Pinterp_rows[i] = row
            pressure = self.min1 + self.delta1 * index_p
        elif self.second_dimension_name == "Pressure":
            if self.Pinterp_rows is None:
                self.Pinterp_rows = np.zeros(self.number1).astype(int)
            for i in range(self.number1):
                row = index_p*self.number1 + i
                self.Pinterp_rows[i] = row
            pressure = self.min2 + self.delta2 * index_p
        else:
            return ValueError("The column of Pressure is not found.")
        
        # pressure and entropy
        for j in range(len(entropies)):
            self.PS_data[PS_rows[0] + j, 0] = entropies[j]
            self.PS_data[PS_rows[0] + j, 1] = pressure

        # extrapolate data 
        entropy_data = self.data[np.ix_(self.Pinterp_rows), col_entropy]

        if debug:
            print("pressure: ")
            print(pressure)
            print("self.Pinterp_rows")
            print(self.Pinterp_rows)
            print("entropy_data: ")
            print(entropy_data) # debug


        for i in range(len(field_names)) :
            field_name = field_names[i]
            col_data = self.header[field_name]['col']
            field_data = self.data[np.ix_(self.Pinterp_rows), col_data]
            if debug:
                print("field_data: ")
                print(field_data)
            temp = np.interp(entropies, entropy_data[0], field_data[0])
            for j in range(len(entropies)):
                self.PS_data[PS_rows[0]+j, i+2] = temp[j]
    
    def InterpolatePressureEntropy(self, entropies, field_names):
        '''
        Interpolate the data to pressure entropy field
        Inputs:
            entropies (list): entropy inputs
            field_names (list): names of field to interpolate
        '''
        # initiate
        n_field = len(field_names)
        n_entropy = len(entropies)
        if self.first_dimension_name == "Pressure":
            n_p = self.number1
        elif self.second_dimension_name == "Pressure":
            n_p = self.number2
        else:
            return ValueError("The column of Pressure is not found.")
        self.PS_data = np.zeros([n_entropy*n_p, n_field + 2]) 

        # call the function to interpolate data for one pressure value
        for index_p in range(n_p):
            PS_rows = [n_entropy * index_p, n_entropy * (index_p + 1)]
            self.InterpolatePressureEntropyByIndex(index_p, entropies, field_names, PS_rows)
        self.PS_number_out1 = n_entropy
        self.PS_number_out2 = n_p


    def PlotHefesto(self):
        '''
        Plot the Hefesto lookup table
    
        Inputs:
            -
        Returns:
            -col_alpagg
        '''
        pass


class CONTROL_FILE():

    def __init__(self):
        '''
        initiation
        Attributes:
            lines: line inputs
            n_raw (int): number of line inputs
            P1, P2, nP: entries for the P / density dimention
            T1, T2, nT: entries for the T / entropy dimention
            useT (int): 0 - use S; 1 - use T
            prm_dir: directory to the parameter files
        '''
        self.lines = []
        self.n_raw = 0
        self.P1 = 0.0
        self.P2 = 0.0
        self.nP = 0
        self.useT = 1
        self.prm_dir = None
        pass

    def ReadFile(self, _path):
        '''
        Inputs:
            _path (str): _path of a control file
        '''
        assert(os.path.isfile(_path))
        with open(_path, 'r') as fin:
            contents = fin.read()
        self.lines = contents.split('\n')
        self.n_raw = len(self.lines)

        # read the first line, dimentions in P, T
        line1 = self.lines[0]
        temp = line1.split(',')
        foo = Utilities.re_neat_word(temp[0])
        self.P1 = float(foo)
        foo = Utilities.re_neat_word(temp[1])
        self.P2 = float(foo)
        foo = Utilities.re_neat_word(temp[2])
        self.nP = int(foo)
        foo = Utilities.re_neat_word(temp[3])
        self.T1 = float(foo)
        foo = Utilities.re_neat_word(temp[4])
        self.T2 = float(foo)
        foo = Utilities.re_neat_word(temp[5])
        self.nT = int(foo)

    def configureP(self, P1, P2, nP):
        '''
        Inputs:
            P1, P2, nP: entries for the P / density dimention
        '''
        self.P1 = P1
        self.P2 = P2
        self.nP = nP
    
    def configureT(self, T1, T2, nT, **kwargs):
        '''
        Inputs:
            T1, T2, nT: entries for the T / entropy dimention
            kwargs:
                use_T: use temperature

        '''
        useT = kwargs.get("useT", 1)
        self.useT = useT
        self.T1 = T1
        self.T2 = T2
        self.nT = nT
    
    def configurePrm(self, prm_dir):
        '''
        Inputs:
            directory of the parameter files
        '''
        self.prm_dir = prm_dir

    def WriteFile(self, _path):
        '''
        Inputs:
            _path (str): _path of a control file to write
        '''
        o_lines = self.lines.copy()
        # first line: P, T dimension
        o_lines[0] = ''
        line1 = self.lines[0]
        temp = line1.split(',')
        foo = "%s,%s,%s" % (self.P1, self.P2, self.nP)
        o_lines[0] += foo
        foo = ",%s,%s,%s" % (self.T1, self.T2, self.nT)
        o_lines[0] += foo
        if self.useT == 1:
            temp1 = 0
        else:
            # use entropy
            temp1 = -2
        foo = ",%s,%s,%s" % (temp1, temp[7], temp[8])
        o_lines[0] += foo
        # path of the parameter files
        if self.prm_dir is not None:
            o_lines[10] = self.prm_dir

        # write file
        with open(_path, 'w') as fout:
            for i in range(self.n_raw):
                line = o_lines[i]
                if i > 0:
                    fout.write('\n')
                fout.write(line)
        assert(os.path.isfile(_path))
        print("Write file %s" % _path)


def ParsePerplexHeader(line):
    '''
    Parse header of a perplex file
    Inputs:
        line (str) - header inputs
    Returns:
        header (list) - a list of header
        unit (list) - a list of unit
    '''
    words = line.split(' ')
    header = []
    unit = []
    for word in words:
        word = Utilities.re_neat_word(word)
        if len(word) > 0:
            temp = word.split(',')
            if len(temp) == 2:
                unit.append(temp[1])
                header.append(temp[0])
            else:
                temp = word.split('(')
                assert(len(temp) == 2)
                unit.append(temp[1].replace(')', ''))
                header.append(temp[0])
    assert(len(header) == len(unit))
    return header, unit


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


def FixFirstDimensionMinor(nddata, number1, **kwargs):
    '''
    fix minor differences in the data in the first dimension
    Inputs:
        nddata: a ndarray of the first dimension data
        number1: number1 i the first dimention
    '''
    limits = kwargs.get('limits', [1e-16, 1e-4])
    for i in range(0, nddata.size):
        if i >= number1:
            if nddata[i-number1] < 1e-32:
                # don't divide a 0 value
                diff = abs(nddata[i] - nddata[i-number1])
            else:
                diff = abs((nddata[i] - nddata[i-number1])/nddata[i-number1])
            if diff > limits[0] and diff < limits[1]:
                nddata[i] = nddata[i-number1]
            elif diff > limits[1]:
                raise ValueError("Two coordinates (%d, %d) in the first dimension have large differences" % (i-number1, i))


def FixSecondDimensionMinor(nddata, number1, **kwargs):
    '''
    fix minor differences in the data in the second dimension
    Inputs:
        nddata: a ndarray of the second dimension data
        number1: number1 in the second dimention
    '''
    limits = kwargs.get('limits', [1e-16, 1e-4])
    for i in range(0, nddata.size):
        if i % number1 != 0:
            if nddata[i-1] < 1e-32:
                # don't divide a 0 value
                diff = abs((nddata[i] - nddata[i-1]))
            else:
                diff = abs((nddata[i] - nddata[i-1])/nddata[i-1])
            if diff > limits[0] and diff < limits[1]:
                nddata[i] = nddata[i-1]
            elif diff > limits[1]:
                raise ValueError("Two coordinates (%d, %d) in the second dimension have large differences" % (i-1, i))


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
    return min, delta, number


def ProcessHefesto(filein, fileout, interval1, interval2):
    # input file
    assert(os.path.isfile(filein))
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadHeFestoTable(filein)
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
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadHeFestoTable(filein)
    LookupTable.Update(first_dimension_name)


def ConvertPS_Table(filein, fileout):
    '''
    Converts a (P, T) -> S lookup table to a (P, S) -> T lookup table
    Inputs:
        filein(str) - path to the input lookup table
        fileout(str) - path to the output lookup table
    '''
    # input file
    assert(os.path.isfile(filein))
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadPerplex(filein, header_rows=4)
    LookupTable.Update()

    # output pressure entropy lookup table
    entropies = np.linspace(550.0, 3300.0, 56)
    field_names = ['Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    output_field_names = ['Entropy', 'Pressure', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.InterpolatePressureEntropy(entropies, field_names)
    LookupTable.OutputPressureEntropyTable(output_field_names, fileout)
    assert(os.path.isfile(fileout))


def PlotBao22_sup_fig1d(ax):
    '''
    Plot data from supplimentary file in Bao 2022, figure 1d
    Inputs:
        ax - a matplotlib axis
    '''
    depths = np.array([268.2926829268293, 291.5989159891599, 298.1029810298103, 308.9430894308943, 322.49322493224935, 337.1273712737127,\
                        340.10840108401084, 341.19241192411926, 355.55555555555554, 371.81571815718155, 395.66395663956644, 402.9810298102981,\
                            406.5040650406504, 409.7560975609756, 410.840108401084, 415.17615176151764, 421.1382113821138, 430.08130081300817,\
                                449.32249322493226, 460.70460704607046, 484.2818428184282, 488.8888888888889, 500.8130081300813, 527.9132791327913,\
                                    535.7723577235772, 542.0054200542006, 556.0975609756097, 599.4579945799458])
    vss = np.array([4.540342298288504, 4.5599022004889935, 4.564792176039116, 4.5721271393643, 4.586797066014666, 4.599022004889972,\
                     4.603911980440094, 4.635696821515889, 4.645476772616133, 4.6601466992665, 4.6797066014669895, 4.68215158924205,\
                        4.83374083129584, 5.002444987775059, 5.061124694376526, 5.066014669926648, 5.073349633251831, 5.085574572127136,\
                            5.107579462102687, 5.119804400977992, 5.134474327628359, 5.1369193154034205, 5.144254278728604, 5.163814180929093,\
                            5.166259168704154, 5.20048899755501, 5.293398533007332, 5.315403422982883])
    ax.plot(depths, vss, "b--", label="Bao 22 sup fig1d")
    ax.set_title("Compare Vs from HeFesto (S = 2.5731454607326154)")
    ax.set_xlabel("Depth (km)")
    ax.set_ylabel("Vs (m/s)")


def CompareHeFestoVs(_path, **kwargs):
    '''
    Compare HeFesto outputs to standard results in paper
    kwargs:
        o_dir: output directory
    '''
    o_dir = kwargs.get("o_dir", RESULT_DIR)
    fig, ax = plt.subplots()
    # plot results from 22 paper
    PlotBao22_sup_fig1d(ax)
    # 
    LookupTable=LOOKUP_TABLE()
    LookupTable.ReadRawFort56(_path)
    depths, VSs = LookupTable.export_fort56_vss()
    ax.plot(depths, VSs, "r", label="Hefesto Outputs")
    # plot settings
    ax.grid()
    ax.set_xlim([250.0, 600.0])
    ax.legend()
    # save figures
    fig_path = os.path.join(o_dir, "HeFesto_compare.png")
    if os.path.isfile(fig_path):
        # remove old files
        os.remove(fig_path)
    fig.savefig(fig_path)
    print("Plot figure %s" % fig_path)


def PlotHeFestoProfile(_path0, **kwargs):
    '''
    Plot the profile from HeFesto. Note that the profile needs to be an adiabat
    Inputs:
        _path0 (str): path of the Hefesto Outputs
        kwargs:
            axT (matplotlib axis): axis to plot the temperatures
            ax_density (matplotlib axis): axis to plot the density
    
    '''
    axT = kwargs.get('axT', None)
    axP = kwargs.get('axP', None)
    ax_density = kwargs.get('ax_density', None)

    # get data 
    LookupTable=LOOKUP_TABLE()
    LookupTable.ReadRawFort56(_path0)
    depths_0, Ts_0 = LookupTable.export_temperature_profile()
    _, densities_0 = LookupTable.export_density_profile()
    _, pressures_0 = LookupTable.export_pressure_profile()
    
    # plot temperature
    if axT is not None:
        axT.plot(Ts_0, depths_0, "--", label="HeFesto T0", color="tab:red")
        axT.set_ylabel("Depth [km]")
        axT.set_xlabel("Temperature [K]")
        axT.legend()
 
    # plot pressure
    if axP is not None:
        axP.plot(pressures_0, depths_0, "--", label="HeFesto P0", color="tab:green")
        axP.set_ylabel("Depth [km]")
        axP.set_xlabel("Pressure [GPa]")
        axP.legend()

    # plot density 
    if ax_density is not None:
        ax_density.plot(densities_0*1000.0, depths_0, "--", label="HeFesto Density", color="tab:blue")
        ax_density.set_ylabel("Depth [km]")
        ax_density.set_xlabel("Density [kg/m^3]")
        ax_density.legend()


    
    return depths_0, densities_0, Ts_0


def ComputeBuoyancy(_path0, _path1, **kwargs):
    '''
    Compute Buoyancy from two Hefesto Outputs
    Inputs
        _path0 (str): path of the first file
        _path1 (str): path of the second file
        kwargs:
            axT (matplotlib axis): axis to plot the temperatures
            ax_density_ratio (matplotlib axis): axis to plot the density differences
            ax_buoy (matplotlib axis): axis to plot the buoyancy ratio
    '''
    n_depth = 1000
    g = 9.8
    axT = kwargs.get('axT', None)
    ax_density_ratio = kwargs.get('ax_density_ratio', None)
    ax_buoy = kwargs.get('ax_buoy', None)
    # read data
    LookupTable=LOOKUP_TABLE()
    LookupTable.ReadRawFort56(_path0)
    depths_0, densities_0 = LookupTable.export_density_profile()
    _, alphas_0 = LookupTable.export_thermal_expansivity_profile()
    _, Ts_0 = LookupTable.export_temperature_profile()
    min_depth0 = depths_0[0]
    max_depth0 = depths_0[-1]
    LookupTable.ReadRawFort56(_path1)
    depths_1, densities_1 = LookupTable.export_density_profile()
    _, alphas_1 = LookupTable.export_thermal_expansivity_profile()
    _, Ts_1 = LookupTable.export_temperature_profile()
    min_depth1 = depths_1[0]
    max_depth1 = depths_1[-1]
    # choose the mutual range in the data
    min_depth = np.max(np.array([min_depth0, min_depth1]))
    max_depth = np.min(np.array([max_depth0, max_depth1]))
    
    # interpolate data
    DensityFunc0 = interpolate.interp1d(depths_0, densities_0, assume_sorted=True, fill_value="extrapolate")
    AlphaFunc0 = interpolate.interp1d(depths_0, alphas_0, assume_sorted=True, fill_value="extrapolate")
    TFunc0 = interpolate.interp1d(depths_0, Ts_0, assume_sorted=True, fill_value="extrapolate")
    DensityFunc1 = interpolate.interp1d(depths_1, densities_1, assume_sorted=True, fill_value="extrapolate")
    AlphaFunc1 = interpolate.interp1d(depths_1, alphas_1, assume_sorted=True, fill_value="extrapolate")
    TFunc1 = interpolate.interp1d(depths_1, Ts_1, assume_sorted=True, fill_value="extrapolate")

    # compute buoyancy and buoyancy number
    # the buoyancy number is computed with buoyancy / density1
    # density1 is chosen instead of density0 to simulate the buoyancy ratio of density0
    depths = np.linspace(min_depth, max_depth, n_depth)
    diff_densities = np.zeros(n_depth)
    buoyancies = np.zeros(n_depth)
    density_ratios = np.zeros(n_depth)
    for i in range(n_depth):
        # get values at depth
        depth = depths[i]
        alpha = AlphaFunc1(depth)
        alpha = np.min(np.array([alpha, 5.0]))
        # print("alpha: ", alpha) # debug
        density0 = DensityFunc0(depth) * 1000.0
        density1 = DensityFunc1(depth) * 1000.0
        T0 = TFunc0(depth)
        T1 = TFunc1(depth)
        diff_density = density1 - density0
        # print("diff_density: ", diff_density) # debug
        diff_densities[i] = diff_density
        buoyancy =  - diff_density * g
        buoyancies[i] = buoyancy
        density_ratios[i] = diff_density / density1

    # plot temperature
    if axT is not None:
        axT.plot(Ts_0, depths_0, "--", label="HeFesto T0", color="tab:red")
        axT.plot(Ts_1, depths_1, "--", label="HeFesto T1", color="tab:blue")
        axT.set_ylabel("Depth [km]")
        axT.set_xlabel("Temperature [K]")
        axT.legend()
    # plot density ratio
    if ax_density_ratio is not None:
        ax_density_ratio.plot(density_ratios, depths, "--", label="HeFesto Density Ratio", color="tab:red")
        ax_density_ratio.set_ylabel("Depth [km]")
        ax_density_ratio.set_xlabel("Density Ratio")
        ax_density_ratio.legend()
    # plot buoyancy
    if ax_buoy is not None:
        ax_buoy.plot(buoyancies, depths, "--", label="HeFesto buoyancy", color="tab:blue")
        ax_buoy.set_ylabel("Depth [km]")
        ax_buoy.set_xlabel("Buoyancy [N/m^3]")
        ax_buoy.legend()

    # return variables
    return depths, buoyancies, density_ratios


def DistributeParallelControl(hefesto_dir, o_dir, case_name, nproc, control_path,\
                              T1, T2, nT, P1, P2, nP, T_variable, slurm_file_base, prm_dir, split_by):
    '''
    Generate controls file for running HeFesto in parallel
    Inputs;
        json_opt (dict or json file): inputs
    '''
    assert(split_by in ["P", "T"])

    if T_variable == "temperature":
        useT = 1
    else:
        useT = 0
    
    # make directory
    case_dir = os.path.join(o_dir, case_name)
    if os.path.isdir(case_dir):
        foo = input("Case directory %s exist, remove? [y/n]" % case_dir)
        if foo == 'y':
            rmtree(case_dir)
        else:
            "Terminating"
            exit(0)
    
    # read file
    ControlFile = CONTROL_FILE()
    ControlFile.ReadFile(control_path)
    ControlFile.configurePrm(prm_dir)

    # P Ranges
    p_ranges = []
    T_ranges = []
    if split_by == "P":
        p_interval = (P2 - P1) / nP
        for i in range(nproc):
            if nP >= nproc and nproc > 1:
                foo = nP // nproc
                foo1 = nP % nproc
                if i == nproc - 1:
                    p_range = [P1 + p_interval * foo*i, P2, int(foo + foo1 + 1)]
                else:
                    p_range = [P1 + p_interval * foo * i, P1 + p_interval * foo * (i+1), int(foo + 1)]
            else:
                p_range = [P1, P2, nP]
            p_ranges.append(p_range)
    elif split_by == "T":
        t_interval = (T2 - T1) / nT
        for i in range(nproc):
            if nT >= nproc and nproc > 1:
                foo = nT  // nproc
                foo1 = nT % nproc
                if i == nproc - 1:
                    t_range = [T1 + t_interval * foo*i, T2, int(foo + foo1 + 1)]
                else:
                    t_range = [T1 + t_interval * foo * i, T1 + t_interval * foo * (i+1), int(foo + 1)]
            else:
                t_range = [T1, T2, nT]
            T_ranges.append(t_range)

    # generate cases 
    os.mkdir(case_dir)
    # make subdirectories
    exe_path = os.path.join(hefesto_dir, "main")
    sh_path = os.path.join(case_dir, "configure.sh")
    sh_contents = "#!/bin/bash\n"
    for iproc in range(nproc):
        sub_dir_name = "sub_%04d" % iproc
        sub_dir = os.path.join(case_dir, sub_dir_name)
        os.mkdir(sub_dir)
        # append new line to sh file 
        o_exe_path = os.path.join(sub_dir_name, "main")
        new_line = "cp %s %s" % (exe_path, o_exe_path)
        sh_contents += (new_line + '\n')
        # configure P, T
        print("split_by: ", split_by) # debug
        if split_by == "P":
            p_range = p_ranges[iproc]
            ControlFile.configureP(p_range[0], p_range[1], p_range[2])
        else:
            ControlFile.configureP(P1, P2, nP)
        if split_by == "T":
            T_range = T_ranges[iproc]
            ControlFile.configureT(T_range[0], T_range[1], T_range[2], useT=useT)
        else:
            ControlFile.configureT(T1, T2, nT, useT=useT)
        # write file
        temp = os.path.join(sub_dir, "control")
        ControlFile.WriteFile(temp)

    # generate bash file 
    with open(sh_path, 'w') as fin:
        fin.write(sh_contents)
    
    # generate the slurm file
    slurm_file_path = os.path.join(case_dir, "job.sh")
    SlurmOperator = ParsePrm.SLURM_OPERATOR(slurm_file_base)
    SlurmOperator.SetAffinity(1, nproc, 1)
    SlurmOperator.ResetCommand()
    SlurmOperator.SetName(case_name)
    SlurmOperator.SetModule(['intel-oneapi-mkl/2022.2.1'], [])

    SlurmOperator.SetTimeByHour(300)
    # generate the command to run
    extra_contents = ""
    temp = "subdirs=("
    for iproc in range(nproc):
        if iproc > 0:
            temp += " "
        sub_dir_name = "sub_%04d" % iproc
        temp += "\"%s\"" % sub_dir_name
    temp += ")\n"
    extra_contents += temp

    temp = """
for subdir in ${subdirs[@]}; do
        cd ${subdir}
        srun --exclusive --ntasks 1 ./main control &
        cd ..
done
wait
"""
    extra_contents += temp
    SlurmOperator.SetExtra(extra_contents)
    SlurmOperator(slurm_file_path)
    print("%s: make new case %s" % (Utilities.func_name(), case_dir))


def CreateHeFESToCaseFromJSON(json_option):
    '''
    Create HeFESTo case from a json file
    Inputs:
        json_option (str or dict): either a file to import or a dict of
            options
    '''
    HeFESTo_Opt = HEFESTO_OPT()
    # read in json options
    if type(json_option) == str:
        if not os.access(json_option, os.R_OK):
            raise FileNotFoundError("%s doesn't exist" % json_option)
        HeFESTo_Opt.read_json(json_option)
    elif type(json_option) == dict:
        HeFESTo_Opt.import_options(json_option)
    else:
        raise TypeError("Type of json_opt must by str or dict")
    
    DistributeParallelControl(*HeFESTo_Opt.to_distribute_parallel_control())

    # copy the json file
    case_dir = HeFESTo_Opt.case_dir()
    json_o_path = os.path.join(case_dir, "case.json")
    if type(json_option) == str:
        copy2(json_option, json_o_path)
    elif type(json_option) == dict:
        with open(json_o_path, 'w') as fout:
            json.dump(fout)
    

def AssembleParallelFiles(case_dir):
    '''
    Inputs:
        case_dir: path to the case directory
    '''
    # screen outputs
    print("start AssembleParallelFiles")

    # load options
    json_file = os.path.join(case_dir, "case.json")
    assert(os.path.isfile(json_file))
    with open(json_file, 'r') as fin:
        case_opt = json.load(fin)
    nproc = case_opt['nproc']
    nP = case_opt['P']['nP']
    # screen outputs
    print("nproc: ", nproc)

    # read file contents 
    # read first file
    sub_dir_name = "sub_0000"
    sub_dir = os.path.join(case_dir, sub_dir_name)
    output_dir = os.path.join(case_dir, "output")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    fort_56_o_path = os.path.join(output_dir, "fort.56")
    # read subsequent files
    fort_56_path = os.path.join(sub_dir, "fort.56")
    with open(fort_56_path, 'r') as fin:
        lines = fin.readlines()
    for i in range(1, nproc):
        sub_dir_name = "sub_%04d" % i
        sub_dir = os.path.join(case_dir, sub_dir_name)
        fort_56_path = os.path.join(sub_dir, "fort.56")
        with open(fort_56_path, 'r') as fin:
            temp = fin.readlines()
        for i in range(nP+1,len(temp)):
            lines.append(temp[i])

    # write file 
    with open(fort_56_o_path, 'w') as fout:
        for _line in lines:
            fout.write(_line)
    
    print("File generated:", fort_56_o_path) # debug


# todo_para
def PlotHeFestoBuoyancy(ifile):
    '''
    Inputs:
        ifile: input file
    '''
    # read input file
    assert(os.path.isfile(ifile))
    with open(ifile, 'r') as fin:
        lines = fin.readlines()
    dir0 = Utilities.re_neat_word(lines[0])
    path0 = os.path.join(dir0, "fort.56")
    dir1 = Utilities.re_neat_word(lines[1])
    path1 = os.path.join(dir1, "fort.56")
    o_dir = Utilities.re_neat_word(lines[2])

    fig = plt.figure(tight_layout=True, figsize=(15, 5))  # plot of wallclock
    gs = gridspec.GridSpec(1, 3)
    axT = fig.add_subplot(gs[0, 0])
    ax_density_ratio = fig.add_subplot(gs[0, 1])
    ax_buoy = fig.add_subplot(gs[0, 2])
    ComputeBuoyancy(path0, path1, axT=axT, ax_density_ratio=ax_density_ratio, ax_buoy=ax_buoy)
    axT.invert_yaxis()
    ax_density_ratio.invert_yaxis()
    ax_buoy.invert_yaxis()
    axT.grid()
    ax_density_ratio.grid()
    ax_buoy.grid()

    o_path = os.path.join(o_dir, "buoyancy.pdf")
    fig.savefig(o_path) # debug
    print("Figure generated: ", o_path)


def convert_mol_fraction(comps):
    '''
    given a composition array of mol%, compute the mol atom %
    Inputs:
        comps - inputs of array of mol%
    Returns:
        comps_atom - inputs of array of atom mol %
    '''
    assert(len(comps) == 6) # SiO2, MgO, FeO, CaO, Al2O3, Na2O
    mol_total = comps[0] + comps[1] + comps[2] + comps[3] + comps[4] + comps[5]
    assert((100.0 - mol_total)/100.0 < 1e-3)
    mol_atom_total = comps[0] + comps[1] + comps[2] + comps[3] + 2 * comps[4] + 2 * comps[5]
    comps_atom = [0 for i in range(6)]
    comps_atom[0] = comps[0] / mol_atom_total
    comps_atom[1] = comps[1] / mol_atom_total
    comps_atom[2] = comps[2] / mol_atom_total
    comps_atom[3] = comps[3] / mol_atom_total
    comps_atom[4] = 2*comps[4] / mol_atom_total
    comps_atom[5] = 2*comps[5] / mol_atom_total
    return comps_atom


def ExchangeDimensions(indexes, number_out1, number_out2):
    '''
    exchange the indexing for the 1st and the 2nd dimensions
    '''
    assert(indexes.ndim==1)
    ixx = np.zeros(indexes.shape, dtype=int)
    i = 0
    for i1 in range(number_out1):
        for j2 in range(number_out2):
            ixx[i] = j2 * number_out1 + i1
            i += 1
    ex_indexes = indexes[np.ix_(ixx)] 

    return ex_indexes
    

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

    elif _commend == "convert_ps":
        # Convert to (P, S) -> T lookup table
        ConvertPS_Table(arg.inputs, arg.outputs)

    elif _commend == "compare_hefesto_outputs":
        # Compare outputs to standard
        CompareHeFestoVs(arg.inputs)

    elif _commend == "create_case":
        CreateHeFESToCaseFromJSON(arg.inputs)

    elif _commend == "assemble_parallel_files":
        AssembleParallelFiles(arg.inputs)

    elif _commend == "plot_hefesto_buoyancy":
        # todo_para 
        PlotHeFestoBuoyancy(arg.inputs)


# run script
if __name__ == '__main__':
    main( )