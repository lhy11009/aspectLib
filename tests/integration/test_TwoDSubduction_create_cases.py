import os
import json
import filecmp
import shilofue.TwoDSubduction0.CreateCases as CreateCases
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
from shutil import rmtree

ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
test_dir = '.test'


def test_generate_case():
    # Test 1: generate case with Lower Mantle viscosity
    parse_dir = os.path.join(test_source_dir, 'TwoDSubduction', 'parse')
    _test_prm_file = os.path.join(parse_dir, 'base.prm')
    _config_file = os.path.join(parse_dir, 'base_config.json')
    _case_dir = os.path.join('.test', 'ULV3.000e+01')  # case name is 'ULV3.000e+01'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(parse_dir, 'standard_base.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParsePrm.ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _config = json.load(fin)
    # Initiate Class MyCase
    MyCase = CreateCases.MYCASE(_inputs, config=_config)
    # call __call__ function
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    # add operations
    parse_operations = CreateCases.MY_PARSE_OPERATIONS()
    MyCase(parse_operations, dirname='.test', extra=_extra)
    # Assertions
    assert(os.path.isfile(_prm_file))
    assert(filecmp.cmp(_standard_prm_file, _prm_file))
    
    # Test 2: generate case with Lower Mantle viscosity with changed mesh_refinement
    _test_prm_file = os.path.join(parse_dir, 'base.prm')
    _config_file = os.path.join(parse_dir, 'base_config_1.json')
    _case_dir = os.path.join('.test', 'ULV1.000e+02testIAR8')  # case name is 'ULV3.000e+01_testIAR8'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(parse_dir, 'standard_base_1.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParsePrm.ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
        _config = _json_inputs['config']
        _test = _json_inputs.get('test', {})
    # Initiate Class MyCase
    MyCase = CreateCases.MYCASE(_inputs, config=_config, test=_test)
    # call __call__ function
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    parse_operations = CreateCases.MY_PARSE_OPERATIONS()
    MyCase(parse_operations, dirname='.test', extra=_extra)
    # Assertions
    assert(os.path.isfile(_prm_file))
    assert(filecmp.cmp(_standard_prm_file, _prm_file))
    
    # Test 3: generate case with non_liear rheology and test solver
    _test_prm_file = os.path.join(parse_dir, 'non_linear_1e18.prm')
    _config_file = os.path.join(parse_dir, 'non_linear_1e18_config.json')
    _case_dir = os.path.join('.test', 'ULV3.000e+01testC4.000e-01MLT9.000e-01NST1.000e-04SBR10')  # case name is 'ULV3.000e+01_testIAR8'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(parse_dir, 'standard_non_linear_1e18.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParsePrm.ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
        _config = _json_inputs['config']
        _test = _json_inputs.get('test', {})
        _extra = _json_inputs.get('extra', {})
    # Initiate Class MyCase
    MyCase = CreateCases.MYCASE(_inputs, config=_config, test=_test, extra=_extra)
    # call __call__ function
    parse_operations = CreateCases.MY_PARSE_OPERATIONS()
    MyCase(parse_operations, dirname='.test', extra=_extra)
    # Assertions
    assert(os.path.isfile(_prm_file))
    assert(filecmp.cmp(_standard_prm_file, _prm_file))
    
    # Test 4: generate case with phase transitions on all compositions and a eclogite transition of crustal layer
    _test_prm_file = os.path.join(parse_dir, 'crust_terminate.prm')
    _config_file = os.path.join(parse_dir, 'crust_terminate_config.json')
    _case_dir = os.path.join('.test', 'crust_terminateULV1.000e+01')  # case name is 'ULV3.000e+01_testIAR8'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(parse_dir, 'crust_terminate_standard.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParsePrm.ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
        _config = _json_inputs['config']
        _test = _json_inputs.get('test', {})
        _extra = _json_inputs.get('extra', {})
    # Initiate Class MyCase
    MyCase = CreateCases.MYCASE(_inputs, config=_config, test=_test, extra=_extra)
    # call __call__ function
    parse_operations = CreateCases.MY_PARSE_OPERATIONS()
    MyCase(parse_operations, basename="crust_terminate", dirname='.test', extra=_extra)
    # Assertions
    assert(os.path.isfile(_prm_file))
    assert(filecmp.cmp(_standard_prm_file, _prm_file))


def test_generate_group():
    # test 1, generate a group
    parse_dir = os.path.join(test_source_dir, 'TwoDSubduction', 'parse')
    _test_prm_file = os.path.join(parse_dir, 'base.prm')
    _config_file = os.path.join(parse_dir, 'base_config_group.json')
    _odir = os.path.join(test_dir, 'test_group')
    if os.path.isdir(_odir):
        # remove older files
        rmtree(_odir)
    os.mkdir(_odir)
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParsePrm.ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
    # Initial test group
    MyGroup = Parse.GROUP_CASE(CreateCases.MYCASE, _inputs, _json_inputs)
    # Call __call__ to generate cases
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    _operations = ['LowerMantle', "MeshRefinement"]  # operations to do
    parse_operations = CreateCases.MY_PARSE_OPERATIONS()
    MyGroup(parse_operations, _odir, operations=_operations, extra=_extra)
    # Assertions
    _case_names = ['ULV1.000e+02testIAR6', 'ULV1.000e+02testIAR8', 'ULV3.000e+01testIAR6','ULV3.000e+01testIAR8']
    for _case_name in _case_names:
        _case_dir = os.path.join(_odir, _case_name)  # case name is 'ULV3.000e+01_testIAR8'
        _prm_file = os.path.join(_case_dir, 'case.prm')
        assert(os.path.isfile(_prm_file))