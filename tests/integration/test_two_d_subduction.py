import os
import json
from shutil import rmtree
from shilofue.TwoDSubduction import MYCASE
from shilofue.Parse import ParseFromDealiiInput, GROUP_CASE

test_source_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
test_dir = '.test'


def test_generate_case():
    # Test 1: generate case with Lower Mantle viscosity
    _test_prm_file = os.path.join(test_source_dir, 'TwoDSubduction', 'twoDSubduction_base.prm')
    _config_file = os.path.join(test_source_dir, 'TwoDSubduction', 'config.json')
    _case_dir = os.path.join('.test', 'ULV3.000e+01')  # case name is 'ULV3.000e+01'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(test_source_dir, 'TwoDSubduction', 'standard_case.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _config = json.load(fin)
    # Initiate Class MyCase
    MyCase = MYCASE(_inputs, config=_config)
    # call __call__ function
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    _operations = ['LowerMantle']  # operations to do
    MyCase(dirname='.test', extra=_extra, operations=_operations)
    # Assertions
    assert(os.path.isfile(_prm_file))
    with open(_standard_prm_file, 'r') as standard_fin:
        # chech the content of the prm file
        with open(_prm_file, 'r') as fin:
            assert(fin.read() == standard_fin.read())
    
    # Test 2: generate case with Lower Mantle viscosity with changed mesh_refinement
    _test_prm_file = os.path.join(test_source_dir, 'TwoDSubduction', 'twoDSubduction_base.prm')
    _config_file = os.path.join(test_source_dir, 'TwoDSubduction', 'config2.json')
    _case_dir = os.path.join('.test', 'ULV1.000e+02testIAR8')  # case name is 'ULV3.000e+01_testIAR8'
    if os.path.isdir(_case_dir):
        # remove older files
        rmtree(_case_dir)
    _prm_file = os.path.join(_case_dir, 'case.prm')
    _standard_prm_file = os.path.join(test_source_dir, 'TwoDSubduction', 'standard_case2.prm')
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
        _config = _json_inputs['config']
        _test = _json_inputs.get('test', {})
    # Initiate Class MyCase
    MyCase = MYCASE(_inputs, config=_config, test=_test)
    # call __call__ function
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    _operations = ['LowerMantle', "MeshRefinement"]  # operations to do
    MyCase(dirname='.test', extra=_extra, operations=_operations)
    # Assertions
    assert(os.path.isfile(_prm_file))
    with open(_standard_prm_file, 'r') as standard_fin:
        # chech the content of the prm file
        with open(_prm_file, 'r') as fin:
            assert(fin.read() == standard_fin.read())


def test_generate_group():
    # test 1, generate a group
    _test_prm_file = os.path.join(test_source_dir, 'TwoDSubduction', 'twoDSubduction_base.prm')
    _config_file = os.path.join(test_source_dir, 'TwoDSubduction', 'config_group.json')
    _odir = os.path.join(test_dir, 'test_group')
    if os.path.isdir(_odir):
        # remove older files
        rmtree(_odir)
    os.mkdir(_odir)
    # Read files
    with open(_test_prm_file, 'r') as fin:
        _inputs = ParseFromDealiiInput(fin)
    with open(_config_file, 'r') as fin:
        _json_inputs = json.load(fin)
    # Initial test group
    MyGroup = GROUP_CASE(MYCASE, _inputs, _json_inputs)
    # Call __call__ to generate cases
    _extra = {'T660': 1663.0, 'P660': 21e9, 'LowerV': 1.5e-6}  # extra configuration
    _operations = ['LowerMantle', "MeshRefinement"]  # operations to do
    MyGroup(_odir, operations=_operations, extra=_extra)
    # Assertions
    _case_names = ['ULV1.000e+02testIAR6', 'ULV1.000e+02testIAR8', 'ULV3.000e+01testIAR6','ULV3.000e+01testIAR8']
    for _case_name in _case_names:
        _case_dir = os.path.join(_odir, _case_name)  # case name is 'ULV3.000e+01_testIAR8'
        _prm_file = os.path.join(_case_dir, 'case.prm')
        assert(os.path.isfile(_prm_file))