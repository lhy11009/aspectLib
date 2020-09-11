import math
import shilofue.Utilities as Utilities


def test_read_header():
    '''
    Test the ReadHeader function
    Asserts:
        key == keywords in header, with ' ' substituded by '_'
        cols == columes in header - 1
        units == units in header or None
    '''
    _texts = ['# 1: Time step', '# 2: Time (years)', '# 12: Iterations for composition solver 1']
    _header = Utilities.ReadHeader(_texts)  # call function
    assert(_header['Time_step']['col'] == 0)
    assert(_header['Time_step']['unit'] == None)
    assert(_header['Time']['col'] == 1)
    assert(_header['Time']['unit'] == 'years')
    assert(_header['Iterations_for_composition_solver_1']['col'] == 11)
    assert(_header['Iterations_for_composition_solver_1']['unit'] == None)


def test_read_header2():
    '''
    Test the ReadHeader2 function
    Asserts:
        key == keywords in header, with ' ' substituded by '_'
        cols == columes in header - 1
        units == units in header or None
    '''
    _texts = [' Time_step   Time    Iterations_for_composition_solver_1']
    _header = Utilities.ReadHeader2(_texts)  # call function
    assert(_header['Time_step']['col'] == 0)
    assert(_header['Time_step']['unit'] == None)
    assert(_header['Time']['col'] == 1)
    assert(_header['Time']['unit'] == None)
    assert(_header['Iterations_for_composition_solver_1']['col'] == 2)
    assert(_header['Iterations_for_composition_solver_1']['unit'] == None)


def test_re_count_indent():
    '''
    Tests the re_count_indent funtion
    Asserts:
        indent = right value
    '''
    # test 1
    _pattern = ' \tfoo'
    _indent = Utilities.re_count_indent(_pattern)
    assert(_indent == 5)
    # test 2
    _pattern = '    foo'
    _indent = Utilities.re_count_indent(_pattern)
    assert(_indent == 4)


def test_ggr2cart2():
    '''
    Tests the ggr2cart2 funtion
    Asserts:
        return value is correct
    '''
    # test 1
    lon = 0.0
    r = 1.0
    x, y = Utilities.ggr2cart2(lon, r)
    # x = 1.0, y = 0.0
    assert(abs(x-1.0)/1.0 < 1e-8)
    assert(y/1.0 < 1e-8)
    # test 2
    lon = math.pi / 4.0
    r = 1.0
    x, y = Utilities.ggr2cart2(lon, r)
    # x = sqrt(2)/2, y = sqrt(2)/2
    temp = 2.0**0.5 / 2.0
    assert(abs(x-temp)/temp < 1e-8)
    assert(abs(y-temp)/temp < 1e-8)


def test_ggr2cart():
    '''
    Tests the ggr2cart funtion
    Asserts:
        return value is correct
    '''
    # test 1
    lon = 0.0
    lat = 0.0
    r = 1.0
    # x = 1.0, y = z = 0.0
    x, y, z = Utilities.ggr2cart(lat, lon, r)
    assert(abs(x-1.0)/1.0 < 1e-8)
    assert(y/1.0 < 1e-8)
    assert(z/1.0 < 1e-8)
    # test 2
    lon = - math.pi / 4.0
    lat = math.pi / 4.0
    r = 1.0
    # x = 0.5 , y = -0.5,  z = sqrt(2)/2
    x, y, z = Utilities.ggr2cart(lat, lon, r)
    assert(abs(x-0.5)/0.5 < 1e-8)
    assert(abs((y+0.5)/0.5) < 1e-8)
    temp = 2.0**0.5 / 2.0
    assert(abs(z-temp)/temp < 1e-8)