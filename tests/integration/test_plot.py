import os
import numpy as np
from shilofue import Plot
from shilofue.Utilities import UNITCONVERT


def test_plot_statistics():
    '''
    A test on ploting statistics results
    '''
    if(os.path.isfile('Statistics.pdf')):
        # remove previous files
        os.remove('Statistics.pdf')
    # test_file = 'fixtures/statistics'
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'statistics')
    assert(os.access(test_file, os.R_OK))
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    Statistics = Plot.STATISTICS_PLOT('Statistics', unit_convert=UnitConvert)
    Statistics(test_file, fileout='./Statistics.pdf')
    assert(os.path.isfile('Statistics.pdf'))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished

def test_plot_depth_average():
    '''
    A test on ploting depth average results
    '''
    if(os.path.isfile('DepthAverage.pdf')):
        # remove previous files
        os.remove('DepthAverage.pdf')
    # test_file = 'fixtures/statistics'
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'depth_average.txt')
    assert(os.access(test_file, os.R_OK))
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    # plot statistics ouput #####
    DepthAverage = Plot.DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert)
    DepthAverage(test_file, fileout='./DepthAverage.pdf')
    assert(DepthAverage.time_step_length == 50)
    assert(DepthAverage.time_step_indexes[-1][-1] == 376)
    assert(abs(DepthAverage.time_step_times[0]-0.0)<1e-6)
    assert(abs(DepthAverage.time_step_times[-1]-2.63571e+06)/2.63571e+06 < 1e-6)
    assert(os.path.isfile('DepthAverage.pdf'))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished
