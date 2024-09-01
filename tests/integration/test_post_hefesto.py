# -*- coding: utf-8 -*-
r"""Test for foo.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_foo.py

descriptions:
    every function is a separate test, combined usage with pytest module
"""


import os
# import pytest
import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
from shilofue.PostHefesto import *
from matplotlib import gridspec
from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'post_hefesto')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_ExchangeDimensions():
    '''
    test function ExchangeDimensions
    '''
    input_file = os.path.join(source_dir, "table_index_test.txt")
    indexes = np.loadtxt(input_file)
    number_out1 = 26
    number_out2 = 4
    ex_indexes = ExchangeDimensions(indexes, number_out1, number_out2)
    # save result to file
    output_file = os.path.join(test_dir, "exchange_dimensions_output.txt")
    if os.path.isfile(output_file):
        os.remove(output_file)
    with open(output_file, 'w') as fout:
        np.savetxt(fout, ex_indexes)
    assert(os.path.isfile(output_file))
    # compare
    output_std_file = os.path.join(source_dir, "exchange_dimensions_output_std.txt")
    assert(os.path.isfile(output_std_file))
    assert(filecmp.cmp(output_file, output_std_file))


def test_process_hefesto_fort56_PS_exchange_dimension():
    '''
    Test processing hefesto table and exchange the 1st and 2nd dimensions
    Asserts:
    '''
    input_file = os.path.join(source_dir, "fort.56.PS")
    assert(os.path.isfile(input_file))  # assert there is an existing Hefesto table
    # output paths
    output_file = os.path.join(test_dir, "hefesto_PS_table_from_fort56_exchange_dimension")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "hefesto_PS_table_from_fort56_exchange_dimension_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadRawFort56(input_file)
    # fields to read in
    field_names = ['Pressure', 'Entropy', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.Process(field_names, output_file, interval1=1, interval2=1, second_dimension="Entropy",\
                        fix_coordinate_minor=True, exchange_dimension=True, file_type="structured")
    # assert something 
    assert(os.path.isfile(output_file))
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))


def test_process_hefesto_fort56_PS():
    '''
    Test processing hefesto table
    Asserts:
    '''
    input_file = os.path.join(source_dir, "fort.56.PS")
    assert(os.path.isfile(input_file))  # assert there is an existing Hefesto table

    # test 1
    # output paths
    output_file = os.path.join(test_dir, "hefesto_PS_table_from_fort56")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "hefesto_PS_table_from_fort56_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadRawFort56(input_file)
    # fields to read in
    field_names = ['Pressure', 'Entropy', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.Process(field_names, output_file, interval1=1, interval2=1, second_dimension="Entropy")
    # assert something 
    assert(os.path.isfile(output_file))
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))
    
    # test 2: add fix minor coordinate
    # output paths
    output_file = os.path.join(test_dir, "hefesto_PS_table_from_fort56_fix_minor")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "hefesto_PS_table_from_fort56_fix_minor_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadRawFort56(input_file)
    # fields to read in
    field_names = ['Pressure', 'Entropy', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.Process(field_names, output_file, interval1=1, interval2=1, second_dimension="Entropy", fix_coordinate_minor=True)
    # assert something 
    assert(os.path.isfile(output_file))
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))


def test_process_hefesto_fort56():
    '''
    Test processing hefesto table
    Asserts:
    '''
    input_file = os.path.join(source_dir, "fort.56.PT")
    assert(os.path.isfile(input_file))  # assert there is an existing Hefesto table
    output_file = os.path.join(test_dir, "hefesto_table_from_fort56")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "hefesto_table_from_fort56_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    
    # call processfunction
    LookupTable = LOOKUP_TABLE()
    LookupTable.ReadRawFort56(input_file)
    # fields to read in
    field_names = ['Pressure', 'Temperature', 'Density', 'Thermal_expansivity', 'Isobaric_heat_capacity', 'VP', 'VS', 'Enthalpy']
    LookupTable.Process(field_names, output_file, interval1=1, interval2=1)
    
    # assert something 
    assert(os.path.isfile(output_file))
    
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))


def test_distribute_parallel_control():
    '''
    assert function DistributeParallelControl
    '''
    case_dir = os.path.join(test_dir, "test_hefesto_parallel")
    # remove older results
    if os.path.isdir(case_dir):
        rmtree(case_dir)

    json_file = os.path.join(source_dir, "test_hefesto.json")
    HeFESTo_Opt = HEFESTO_OPT()
    # read in json options
    if type(json_file) == str:
        if not os.access(json_file, os.R_OK):
            raise FileNotFoundError("%s doesn't exist" % json_file)
        HeFESTo_Opt.read_json(json_file)
    elif type(json_file) == dict:
        HeFESTo_Opt.import_options(json_file)
    else:
        raise TypeError("Type of json_opt must by str or dict")
    
    DistributeParallelControl(*HeFESTo_Opt.to_distribute_parallel_control())

    # check directories
    assert(os.path.isdir(case_dir))
    sub0_dir = os.path.join(case_dir, "sub_0000")
    assert(os.path.isdir(sub0_dir))
    control0_path = os.path.join(sub0_dir, "control")
    assert(os.path.isfile(control0_path))
    sub1_dir = os.path.join(case_dir, "sub_0001")
    assert(os.path.isdir(sub1_dir))
    sub2_dir = os.path.join(case_dir, "sub_0002")
    assert(os.path.isdir(sub2_dir))

    # check the control file
    control2_std_path = os.path.join(source_dir, "control_std")
    control2_path = os.path.join(sub2_dir, "control")
    assert(os.path.isfile(control2_path))
    assert(filecmp.cmp(control2_path, control2_std_path))



def test_control_file():
    '''
    test class CONTROL_FILE
    '''
    ifile = os.path.join(source_dir, "control")
    ofile = os.path.join(test_dir, "test_control")
    ControlFile = CONTROL_FILE()
    # read file
    ControlFile.ReadFile(ifile)
    # configure P, T
    ControlFile.configureP(0.0, 30.0, 31)
    # write file
    ControlFile.WriteFile(ofile)
    ofile_std = os.path.join(source_dir, "test_control_std")
    assert(filecmp.cmp(ofile, ofile_std))

    pass


def test_process_hefesto_table():
    '''
    Test processing hefesto table
    Asserts:
    '''
    input_file = os.path.join(source_dir, "example_table")
    assert(os.path.isfile(input_file))  # assert there is an existing Hefesto table
    output_file = os.path.join(test_dir, "example_table_output")
    if (os.path.isfile(output_file)):  # remove old files
        os.remove(output_file)
    output_file_std = os.path.join(source_dir, "example_table_output_std")
    assert(os.path.isfile(output_file_std))  # assert there is an existing standard file
    ProcessHefesto(input_file, output_file, 1, 1)  # process Hefesto output
    assert(os.path.isfile(output_file))  # assert there is output.
    # filecmp
    assert(filecmp.cmp(output_file, output_file_std))

    # assert something
    assert(True)


def test_compute_buoyancy():
    '''
    Test the ComputeBuoyancy function
    Asserts:
        plot utilities work
        values of the buoyancy calculated
    '''
    fort56_file0 = os.path.join(source_dir, "fort.56.0")
    fort56_file1 = os.path.join(source_dir, "fort.56.1")
    fig_path = os.path.join(test_dir, "buoyancy_hefesto.png")
    fig = plt.figure(tight_layout=True, figsize=(6, 10))
    gs = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    if os.path.isfile(fig_path):
        os.remove(fig_path)

    # plot figures
    _, buoyancies, buoyancy_ratios = ComputeBuoyancy(fort56_file0, fort56_file1, odir=test_dir, ax1=ax1, ax2=ax2)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    # save figures
    fig_path = os.path.join(fig_path)
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))

    # check on the buoyancy calculated
    buoyancy440 = -1.1861446688010937
    buoyancy_ratio440 = -0.24262373534926984
    assert(abs(buoyancies[440] - buoyancy440)/buoyancy440 < 1e-6)
    assert(abs(buoyancy_ratios[440] - buoyancy_ratio440)/buoyancy_ratio440 < 1e-6)


# notes

# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

