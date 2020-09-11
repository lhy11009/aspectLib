import pytest
import numpy as np
import shilofue.Utilities as Utilities
from shilofue.Utilities import UNITCONVERT
from importlib import resources

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