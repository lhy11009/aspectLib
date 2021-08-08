import os
import json
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
from shutil import rmtree

_test_dir = ".test"
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(_test_dir):
    # check we have the directory to store test result
    os.mkdir(_test_dir)


def test_change_disc_values():
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_test.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = ParsePrm.ParseFromDealiiInput(fin)
    _config = {'names': [['End time'], ['Material model', 'Visco Plastic', 'Reset corner viscosity constant']],
               'values': ['80.0e6', '1e21']}
    Parse.ChangeDiscValues(inputs, _config['names'], _config['values'])  # call function
    assert(inputs['End time'] == '80.0e6')
    assert(inputs['Material model']['Visco Plastic']['Reset corner viscosity constant'] == '1e21')
    pass

def test_parse_to_new_case():
    # test usage1: 
    #   step a: read a parameter file
    #   step b: initiate a CASE class with parameters
    #   step c: call calss function UpdatePrmDict explicitly to update parameters
    #   step d: call class __class__ function by method='manual' and filename=some_file
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_test.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = ParsePrm.ParseFromDealiiInput(fin)
    Case = Parse.CASE(inputs)
    _names = [['End time'], ['Material model', 'Visco Plastic', 'Reset corner viscosity constant']]
    _values = ['80.0e6', '1e21']
    Case.UpdatePrmDict(_names, _values)
    assert(Case.idict['End time'] == '80.0e6')
    assert(Case.idict['Material model']['Visco Plastic']['Reset corner viscosity constant'] == '1e21')
    # call CASE: __call__ function to generate a output
    _ofile = os.path.join(_test_dir, 'foo.prm')
    if os.path.isfile(_ofile):
        # remove older file
        os.remove(_ofile)
    parse_operations = Parse.PARSE_OPERATIONS()
    Case(parse_operations, method='manual', filename=_ofile)
    assert(os.path.isfile(_ofile))
    # test usage1: 
    #   step a: read a parameter file
    #   step b: initiate a CASE class with parameters and a config dictionary
    #   step d: call class __class__ function by method='auto'
    _odir = os.path.join('.test', 'test_case_by_auto')
    _ofile = os.path.join(_odir, 'case.prm')
    if os.path.isdir(_odir):
        # remove older directory
        rmtree(_odir)
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_test.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = ParsePrm.ParseFromDealiiInput(fin)
    Case = Parse.CASE(inputs, config={})
    parse_operations = Parse.PARSE_OPERATIONS()
    Case(parse_operations, dirname='.test', basename='test_case_by_auto')
    assert(os.path.isfile(_ofile))


def test_get_sub_cases():
    case_dirs = Parse.GetSubCases(test_source_dir)
    assert(case_dirs == 
           ['/home/lochy/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/parse/foo1',
           '/home/lochy/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/parse/foo'])