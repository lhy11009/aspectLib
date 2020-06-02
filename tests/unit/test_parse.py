import pytest
from shilofue.Parse import GetGroupCaseFromDict
from shilofue.Parse import ExpandNamesParameters

def test_build_group_from_dict():
    # test 0-1, raise error because input type is not dict
    with pytest.raises(TypeError) as excinfo:
        _names, _parameters = GetGroupCaseFromDict('foo')
    assert('dictionary' in str(excinfo.value))
    # test 0-2, raise error because some entry is not int, float or str
    _idict = {}
    _idict['test parameter'] = None
    with pytest.raises(TypeError) as excinfo:
        _names, _parameters = GetGroupCaseFromDict(_idict)
    assert('int, float or str' in str(excinfo.value))
    # test 0-3 
    # todo error in ExpandNamesParameters function

    # test 1, corresponds to
    # 'set Test parameter = 0 or 1 or 2'
    # in a .prm file
    _idict={}
    _idict['test parameter'] = [0, 1, 2]
    _names, _parameters = GetGroupCaseFromDict(_idict)
    assert(_names == [['test parameter']])
    assert(_parameters == [[0, 1, 2]])
    _cases_config = ExpandNamesParameters(_names, _parameters)
    # _cases_config contains 3 cases, which is a list of 3
    # each case has a dictionary with two keys: 'names', 'values'
    assert(len(_cases_config) == 3)
    assert(_cases_config[0]['names'] == [['test parameter']])
    assert(_cases_config[1]['names'] == [['test parameter']])
    assert(_cases_config[2]['names'] == [['test parameter']])
    assert(_cases_config[0]['values'] == [0])
    assert(_cases_config[1]['values'] == [1])
    assert(_cases_config[2]['values'] == [2])
    # test 2, correspond to
    #   'Subsection Material model
    #       subsection Visco Plastic
    #           set Reference viscosity = 1e20 or 1e21 or 1e22
    #           set Maximum viscosity = 1e24
    # 
    #   '
    # in a .prm file
    _idict={}
    _idict['Start time'] = [0, 1e6]
    _idict['Material model'] = {'Visco Plastic': {'Reference viscosity': [1e20, 1e21, 1e22], 'Maximum viscosity': [1e24, 1e25]}}
    _names, _parameters = GetGroupCaseFromDict(_idict)
    assert(_names == [['Material model', 'Visco Plastic', 'Maximum viscosity'], ['Material model', 'Visco Plastic', 'Reference viscosity'] ,['Start time']])
    assert(_parameters == [[1e24, 1e25], [1e20, 1e21, 1e22], [0, 1e6]])
    _cases_config = ExpandNamesParameters(_names, _parameters)
    assert(len(_cases_config) == 12)
    for i in range(12):
        assert(_cases_config[i]['names'] == _names)
    assert (_cases_config[0]['values'] == [1e24, 1e20, 0])
    assert (_cases_config[1]['values'] == [1e24, 1e20, 1e6])
    assert (_cases_config[11]['values'] == [1e25, 1e22, 1e6])

