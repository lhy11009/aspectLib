import pytest
import numpy as np
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
    lon = np.pi / 4.0
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
    lon = - np.pi / 4.0
    lat = np.pi / 4.0
    r = 1.0
    # x = 0.5 , y = -0.5,  z = sqrt(2)/2
    x, y, z = Utilities.ggr2cart(lat, lon, r)
    assert(abs(x-0.5)/0.5 < 1e-8)
    assert(abs((y+0.5)/0.5) < 1e-8)
    temp = 2.0**0.5 / 2.0
    assert(abs(z-temp)/temp < 1e-8)


def test_cart2sph():
    '''
    Tests the cart2sph funtion
    Asserts:
        return value is correct
    '''
    # test1
    x, y, z = 0., 1., 0.
    r0, th0, ph0 = 1., np.pi/2.0, np.pi/2.0
    r, th, ph = Utilities.cart2sph(x, y, z)
    assert(abs((r - r0)/r0) < 1e-8)
    assert(abs((th-th0)/th0) < 1e-8)
    assert(abs((ph-ph0)/ph0) < 1e-8)

    # test2
    x, y, z = 3.0**0.5/2.0, 0.5, 1.0 
    r0, th0, ph0 = 2.0**0.5, np.pi/4.0, np.pi/6.0
    r, th, ph = Utilities.cart2sph(x, y, z)
    assert(abs((r - r0)/r0) < 1e-8)
    assert(abs((th-th0)/th0) < 1e-8)
    assert(abs((ph-ph0)/ph0) < 1e-8)


def test_cart2sph2():
    '''
    Tests the cart2sph2 funtion
    Asserts:
        return value is correct
    '''
    # test1
    x, y = 0., 1.
    r0, ph0 = 1., np.pi/2.0
    r, ph = Utilities.cart2sph2(x, y)
    assert(abs((r - r0)/r0) < 1e-8)
    assert(abs((ph-ph0)/ph0) < 1e-8)

    # test2
    x, y = -0.5, 3.0**0.5/2.
    r0, ph0 = 1., 2.0*np.pi/3.0
    r, ph = Utilities.cart2sph2(x, y)
    assert(abs((r - r0)/r0) < 1e-8)
    assert(abs((ph-ph0)/ph0) < 1e-8)


def test_make_2d_array():
    """
    Tests the Make2dArray function
    Asserts:
        return value is of the correct type
    """
    # test 1
    x = 1.0
    y = Utilities.Make2dArray(x)
    assert(y.shape == (1, 1))

    # test 2
    x = [0.0, 1.0]
    y = Utilities.Make2dArray(x)
    assert(y.shape == (1, 2))
    
    # test 3, we flip test 2
    x = [0.0, 1.0]
    y = np.transpose(Utilities.Make2dArray(x))
    assert(y.shape == (2, 1))

    # test 4, raises error if type is incorrect
    with pytest.raises(TypeError) as excinfo:
        x = {'foo': 0.0}
        y = Utilities.Make2dArray(x)
    assert ('Make2dArray, x must be int, float, list or 1d ndarray' in str(excinfo.value))