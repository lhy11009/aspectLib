import pytest
import os
import numpy as np
from shilofue.Utilities import UNITCONVERT
from importlib import resources
from matplotlib import pyplot as plt
from matplotlib import ticker, cm
import shilofue.Utilities as Utilities

_test_dir = '.test'

def test_unit_convert_json():
    '''
    Test UNITCONVERT class from shilofue.Utilities with
    the json file in shilofue.json 'ConvertUnit.json'
    Assert:
        Units are converted correctly;
        Alias could be used;
        a wrong name of unit is handled;
        a wrong convertion is handled 
    '''
    # assert that a wrong file name wil raise error
    with pytest.raises(FileNotFoundError) as excinfo:
        UnitConvert = Utilities.UNITCONVERT(filename='foo')
    assert('foo' in str(excinfo.value))
    # assert of converting units
    UnitConvert = Utilities.UNITCONVERT()
    assert(np.abs((UnitConvert('km', 'm') - 1000.0)/1000.0)<1e-6)
    assert(np.abs((UnitConvert('m', 'km') - 1/1000.0)/(1/1000.0))<1e-6)
    assert(np.abs((UnitConvert('m/year', 'm/s') - 3.1688739e-8)/3.1688739e-8)<1e-6)
    pytest.raises(KeyError, UnitConvert, 'foo', 'm')  # assert that a wrong name of unit is handled
    pytest.raises(KeyError, UnitConvert, 'm', 'foo')  # assert that a wrong name of unit is handled
    pytest.raises(AssertionError, UnitConvert, 'yr', 'km')  # assert that a wrong convertion is handled 
    # assert that missing alias also causes errors
    UnitConvert = Utilities.UNITCONVERT()
    UnitConvert.alias['foo1'] = 'foo'
    with pytest.raises(KeyError) as excinfo:
        UnitConvert('m', 'foo1')
    assert('alias' in str(excinfo.value))
    with pytest.raises(KeyError) as excinfo:
        UnitConvert('foo1', 'm')
    assert('alias' in str(excinfo.value))


def test_average_phase_function_inputs():
    '''
    todo
    Test AveragePhaseFunction
    '''
    # create data array
    nx = 100
    ny = 200
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xx, yy = np.meshgrid(x, y)
    # zz = np.zeros(xx.shape)
    # derive phase_function_value
    average = Utilities.AveragePhaseFunctionInputs(xx, yy)
    zz = Utilities.PhaseFunction(average)
    # plot
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[0].pcolormesh(xx, yy, average)
    ax[1].pcolormesh(xx, yy, zz)
    filename = os.path.join(_test_dir, "AveragePhaseFunction.png")
    fig.tight_layout()
    fig.savefig(filename)

    assert(os.path.isfile(filename))
    pass