import os
import pytest
import json
import filecmp
import numpy as np
import shilofue.Plot as Plot
from shilofue.Utilities import UNITCONVERT
from matplotlib import pyplot as plt


_test_dir = '.test'
if not os.path.isdir(_test_dir):
    os.mkdir(_test_dir)
_test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test-plot')
assert(os.path.isdir(_test_source_dir))


def test_plot_newton_solver():
    '''
    A test on ploting newton solver results
    '''
    test_file = os.path.join(_test_source_dir, 'newton_solver')
    assert(os.access(test_file, os.R_OK))
    
    # plot stepwise ouput
    _ofile_route = os.path.join(_test_dir, 'NewtonSolverStep.pdf')
    _ofile = os.path.join(_test_dir, 'NewtonSolverStep_s0000000.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)

    # use a json file
    json_file = os.path.join(_test_source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('NewtonSolverStep', {})

    NewtonSolverStep = Plot.NEWTON_SOLVER_PLOT('NewtonSolverStep', options=plot_options)
    
    # plot step0
    NewtonSolverStep.GetStep(0)
    NewtonSolverStep(test_file, fileout=_ofile_route)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully
    
    # plot for all steps
    _ofile = os.path.join(_test_dir, 'NewtonSolver.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)

    # use a json file
    json_file = os.path.join(_test_source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('NewtonSolver', {})
    
    NewtonSolver = Plot.NEWTON_SOLVER_PLOT('NewtonSolver', options=plot_options)
    NewtonSolver(test_file, fileout=_ofile)
    # assert that the file is generated successfull
    assert(os.path.isfile(_ofile))


def test_plot_machine_time():
    '''
    A test on ploting newton solver results
    '''
    test_file = os.path.join(_test_source_dir, 'machine_time')
    assert(os.access(test_file, os.R_OK))

    # use a json file
    json_file = os.path.join(_test_source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('MachineTime', {})
    
    # initiate 
    MachineTime = Plot.MACHINE_TIME_PLOT('MachineTime', options=plot_options)

    # get machine time at one step
    step = 35
    machine_time_at_step, number_of_core = MachineTime.GetStepMT(test_file, step)
    machine_time_at_step_std = 256.0
    number_of_core_std = 128.0
    assert(abs(number_of_core - number_of_core_std) < 1.0)
    assert(abs(machine_time_at_step - machine_time_at_step_std) < 1.0)
    
    step = 100
    machine_time_at_step, number_of_core = MachineTime.GetStepMT(test_file, step)
    machine_time_at_step_std = 675.56
    number_of_core_std = 128.0
    assert(abs(number_of_core - number_of_core_std) < 1.0)
    assert(abs(machine_time_at_step - machine_time_at_step_std) < 1.0)

    # plot 
    _ofile = os.path.join(_test_dir, 'MachineTime.pdf')
    MachineTime(test_file, fileout=_ofile)
    # assert that the file is generated successfull
    assert(os.path.isfile(_ofile))


def test_export_file():
    # export
    '''
    A test on ploting newton solver results
    '''
    test_file = os.path.join(_test_source_dir, 'statistics')
    assert(os.access(test_file, os.R_OK))
    _ofile = os.path.join(_test_dir, 'statistics_export')
    _ofile1 = os.path.join(_test_dir, 'statistics_export1')  # partially output rows
    std_file = os.path.join(_test_source_dir, 'statistics_export')
    std_file1 = os.path.join(_test_source_dir, 'statistics_export1')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    if(os.path.isfile(_ofile1)):
        # remove previous files
        os.remove(_ofile1)
    # sort header
    Statistics = Plot.MACHINE_TIME_PLOT('Statistics')
    Statistics.ReadHeader(test_file)
    Statistics.ReadData(test_file)
    cols, names, units = Statistics.SortHeader()
    assert(names[0] == "Time_step_number")
    assert(names[1] == "Time")
    # export
    Statistics.export(_ofile, ["Iterations_for_temperature_solver", "Iterations_for_composition_solver_1"])
    assert(filecmp.cmp(_ofile, std_file))
    # export1
    Statistics.export(_ofile1, ["Iterations_for_temperature_solver", "Iterations_for_composition_solver_1"], rows=[i for i in range(10)])
    assert(filecmp.cmp(_ofile1, std_file1))