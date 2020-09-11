import os
import shilofue.Parse as Parse
from shutil import rmtree

_test_dir = ".test"
if not os.path.isdir(_test_dir):
    # check we have the directory to store test result
    os.mkdir(_test_dir)

def test_parse_from_file():
    # test_file = 'fixtures/parse_test.prm'
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_test.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = Parse.ParseFromDealiiInput(fin)
    assert(inputs['Dimension'] == '2')
    assert(inputs['Use years in output instead of seconds'] == 'true')
    assert(inputs['End time'] == '40.0e6')
    assert(inputs['Additional shared libraries']
           == '/home/lochy/ASPECT_PROJECT/aspect_plugins/subduction_temperature2d/libsubduction_temperature2d.so, /home/lochy/ASPECT_PROJECT/aspect_plugins/prescribe_field/libprescribed_temperature.so')

def test_change_disc_values():
    test_file = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse_test.prm')
    assert(os.access(test_file, os.R_OK))
    with open(test_file, 'r') as fin:
        inputs = Parse.ParseFromDealiiInput(fin)
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
        inputs = Parse.ParseFromDealiiInput(fin)
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
        inputs = Parse.ParseFromDealiiInput(fin)
    Case = Parse.CASE(inputs, config={})
    parse_operations = Parse.PARSE_OPERATIONS()
    Case(parse_operations, dirname='.test', basename='test_case_by_auto')
    assert(os.path.isfile(_ofile))