import os
import pytest
import json
import numpy as np
import shilofue.Plot as Plot
from shilofue.Utilities import UNITCONVERT
from matplotlib import pyplot as plt


_test_dir = '.test'
if not os.path.isdir(_test_dir):
    os.mkdir(_test_dir)
_test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test-plot')
assert(os.path.isdir(_test_source_dir))


def test_plot_statistics():
    '''
    A test on ploting statistics results
    '''
    _ofile = os.path.join(_test_dir, 'Statistics.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    test_file = os.path.join(_test_source_dir, 'statistics')
    assert(os.access(test_file, os.R_OK))

    # use a json file
    json_file = os.path.join(_test_source_dir, 'linear_plot.json')
    assert(os.access(json_file, os.R_OK))
    with open(json_file, 'r') as fin:
        json_options = json.load(fin)
    plot_options = json_options.get('Statistics', {})

    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    Statistics = Plot.STATISTICS_PLOT('Statistics', unit_convert=UnitConvert, options=plot_options)
    Statistics(test_file, fileout=_ofile)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished


    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully


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