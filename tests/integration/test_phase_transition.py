# -*- coding: utf-8 -*-
r"""Test for foo.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_foo.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
# import pytest
# import filecmp  # for compare file contents
# import numpy as np
from shilofue.PhaseTransition import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'phase_transitions')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_parse_phase_input():
    '''
    test pares_phase_input
    Asserts:
        the output string is the same with the standard one
    '''
    file_path = os.path.join(source_dir, 'phases.json')
    assert(os.access(file_path, os.R_OK))
    outputs = ParsePhaseTransitionFile(file_path)
    assert(outputs==\
    "density = background: 3300.0|3394.4|3442.1|3453.2|3617.6|3691.5|3774.7|3929.1, spharz: 3235.0|3372.3|3441.7|3441.7|3680.8|3717.8|3759.4|3836.6, spcrust: 3000.0|3540.0|3613.0|3871.7"\
    )

    # assert something 
    assert(True)

def test_get_entropy_change():
    '''
    test Get_entropy_change
    assert:
        entropy change on phases for pyrolite
    '''
    file_path = os.path.join(source_dir, 'phases.json')
    assert(os.access(file_path, os.R_OK))
    cdpt_opt = CDPT_OPT()
    cdpt_opt.read_json(file_path)
    phase_opt = cdpt_opt.get_compositions()[0]
    lh = Get_entropy_change(*phase_opt.to_get_latent_heat_contribution())
    assert(str(lh)==\
    '[-33.69606859576077, -16.749786990544447, -3.720269594206477, 26.305683625468298, -22.139111832177445, 18.51104889311546, -13.531580172935637]'\
    )

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

