import pytest
import shilofue.Parse as Parse


def test_build_group_from_dict():
    # test 0-1, raise error because input type is not dict
    with pytest.raises(TypeError) as excinfo:
        _names, _parameters = Parse.GetGroupCaseFromDict('foo')
    assert('dictionary' in str(excinfo.value))
    # test 0-2, raise error because some entry is not int, float or str
    _idict = {}
    _idict['test parameter'] = None
    with pytest.raises(TypeError) as excinfo:
        _names, _parameters = Parse.GetGroupCaseFromDict(_idict)
    assert('int, float or str' in str(excinfo.value))
    # test 0-3 
    # todo error in ExpandNamesParameters function

    # test 1, corresponds to
    # 'set Test parameter = 0 or 1 or 2'
    # in a .prm file
    _idict={}
    _idict['test parameter'] = [0, 1, 2]
    _names, _parameters = Parse.GetGroupCaseFromDict(_idict)
    assert(_names == [['test parameter']])
    assert(_parameters == [[0, 1, 2]])
    _cases_config = Parse.ExpandNamesParameters(_names, _parameters)
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
    _names, _parameters = Parse.GetGroupCaseFromDict(_idict)
    assert(_names == [['Material model', 'Visco Plastic', 'Maximum viscosity'], ['Material model', 'Visco Plastic', 'Reference viscosity'] ,['Start time']])
    assert(_parameters == [[1e24, 1e25], [1e20, 1e21, 1e22], [0, 1e6]])
    _cases_config = Parse.ExpandNamesParameters(_names, _parameters)
    assert(len(_cases_config) == 12)
    for i in range(12):
        assert(_cases_config[i]['names'] == _names)
    assert (_cases_config[0]['values'] == [1e24, 1e20, 0])
    assert (_cases_config[1]['values'] == [1e24, 1e20, 1e6])
    assert (_cases_config[11]['values'] == [1e25, 1e22, 1e6])


def test_get_group_case_from_dict1():
    '''
    Test the function GetGroupCaseFromDict1
    '''
    # test1
    # test dictionary
    _idict = {'config': {'foo': [0, 1, 2]},
              'test': {'foo1': [1, 2, 3]}}
    # call function
    _config_tests = Parse.GetGroupCaseFromDict1(_idict)
    # assertion
    _standard_config_tests = [
                              {'config': {'foo': 0}, 'test': {'foo1': 1}},
                              {'config': {'foo': 1}, 'test': {'foo1': 1}},
                              {'config': {'foo': 2}, 'test': {'foo1': 1}},
                              {'config': {'foo': 0}, 'test': {'foo1': 2}},
                              {'config': {'foo': 1}, 'test': {'foo1': 2}},
                              {'config': {'foo': 2}, 'test': {'foo1': 2}},
                              {'config': {'foo': 0}, 'test': {'foo1': 3}},
                              {'config': {'foo': 1}, 'test': {'foo1': 3}},
                              {'config': {'foo': 2}, 'test': {'foo1': 3}}
                             ]
    assert(_config_tests == _standard_config_tests)
    # test2
    # test dictionary
    _idict = {'config': {'foo': [0, 1], 'foo1': [0, 1]},
              'test': {'foo2': [0, 1]}}
    # call function
    _config_tests = Parse.GetGroupCaseFromDict1(_idict)
    # assertion
    _standard_config_tests = [
                              {'config': {'foo': 0, 'foo1': 0}, 'test': {'foo2': 0}},
                              {'config': {'foo': 1, 'foo1': 0}, 'test': {'foo2': 0}},
                              {'config': {'foo': 0, 'foo1': 1}, 'test': {'foo2': 0}},
                              {'config': {'foo': 1, 'foo1': 1}, 'test': {'foo2': 0}},
                              {'config': {'foo': 0, 'foo1': 0}, 'test': {'foo2': 1}},
                              {'config': {'foo': 1, 'foo1': 0}, 'test': {'foo2': 1}},
                              {'config': {'foo': 0, 'foo1': 1}, 'test': {'foo2': 1}},
                              {'config': {'foo': 1, 'foo1': 1}, 'test': {'foo2': 1}}
    ]
    assert(_config_tests == _standard_config_tests)

    # test 3: include sub_group options
    # _idict = {'config': {'foo': [0, 1], 'foo1': [0, 1], 'sub_group1': {['foo2', 'foo3']: [[2, 3], [4, 5]]}},
    #          'test': {'foo2': [0, 1]}}


def test_change_disc_values():
    '''
    Test the function ChangeDiscValues()
    '''
    # test 1
    # initiate the dictionary
    _idict={}
    _idict['Start time'] = 0
    _idict['Material model'] = {'Visco Plastic': {'Reference viscosity': 1e20, 'Maximum viscosity': 1e24}}
    # initiate two lists
    _names = [['Material model', 'Visco Plastic', 'Maximum viscosity'], ['Material model', 'Visco Plastic', 'Reference viscosity'] ,['Start time']]
    _values = [1e25, 1e22, 1e6]
    Parse.ChangeDiscValues(_idict, _names, _values)  # change value within _idict
    assert(_idict['Start time'] == 1e6)
    assert(_idict['Material model'] == {'Visco Plastic': {'Reference viscosity': 1e22, 'Maximum viscosity': 1e25}})


def test_pattern_from_str():
    '''
    Test the function PatterFromStr
    '''
    assert(Parse.PatternFromStr('foo_foo') == 'FF')
    assert(Parse.PatternFromStr('upper_lower_viscosity') == 'ULV')
    pass


def test_pattern_from_value():
    '''
    Test the function PatterFromValue
    '''
    assert(Parse.PatternFromValue('spherical') == 's')
    assert(Parse.PatternFromValue(0.1) == '1.000e-01')
    assert(Parse.PatternFromValue(100) == '100')
    pass
