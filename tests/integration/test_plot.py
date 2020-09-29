import os
import pytest
import numpy as np
import shilofue.Plot as Plot
from shilofue.Utilities import UNITCONVERT


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
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    json_dir = os.path.join(_test_source_dir, 'json')
    Statistics = Plot.STATISTICS_PLOT('Statistics', unit_convert=UnitConvert, json_dir=json_dir)
    Statistics(test_file, fileout=_ofile)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished


def test_plot_depth_average():
    '''
    A test on ploting depth average results
    '''
    _ofile_route = os.path.join(_test_dir, 'DepthAverage.pdf')
    _ofile = os.path.join(_test_dir, 'DepthAverage_t0.00000000e+00.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    # test_file = 'fixtures/statistics'
    test_file = os.path.join(_test_source_dir, 'depth_average.txt')
    assert(os.access(test_file, os.R_OK))
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    DepthAverage = Plot.DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert)
    # test error handling of key word time
    with pytest.raises(TypeError) as _excinfo:
        DepthAverage(test_file, fileout=_ofile_route, time='foo')
    assert(r'type of time' in str(_excinfo.value))
    # similar, but this time error is on an entry
    with pytest.raises(TypeError) as _excinfo:
        DepthAverage(test_file, fileout=_ofile_route, time=['foo'])
    assert(r'type of values in time' in str(_excinfo.value))
    # generate a successful output with time = 0.0
    DepthAverage(test_file, fileout=_ofile_route, time=0.0)
    assert(DepthAverage.time_step_length == 50)
    assert(DepthAverage.time_step_indexes[-1][-1] == 376)
    assert(abs(DepthAverage.time_step_times[0]-0.0) < 1e-6)
    assert(abs(DepthAverage.time_step_times[-1]-2.63571e+06)/2.63571e+06 < 1e-6)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully


def test_plot_newton_solver():
    '''
    todo
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
    json_dir = os.path.join(_test_source_dir, 'json')
    NewtonSolverStep = Plot.NEWTON_SOLVER_PLOT('NewtonSolverStep', json_dir=json_dir)
    # plot step0
    NewtonSolverStep.GetStep(0)
    NewtonSolverStep(test_file, fileout=_ofile_route)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully
    
    # plot for all steps
    _ofile = os.path.join(_test_dir, 'NewtonSolver.pdf')
    if(os.path.isfile(_ofile)):
        # remove previous files
        os.remove(_ofile)
    json_dir = os.path.join(_test_source_dir, 'json')
    NewtonSolver = Plot.NEWTON_SOLVER_PLOT('NewtonSolver', json_dir=json_dir)
    NewtonSolver(test_file, fileout=_ofile)
    assert(os.path.isfile(_ofile))  # assert that the file is generated successfully