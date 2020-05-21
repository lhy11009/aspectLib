import os
import numpy as np
from shilofue import Plot


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
    # plot statistics ouput #####
    Statistics = Plot.STATISTICS()
    Statistics(test_file)
    assert(os.path.isfile('Statistics.pdf'))  # assert that the file is generated successfully
    # os.remove('Statistics.pdf')  # remove this file after finished
