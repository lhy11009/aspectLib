# -*- coding: utf-8 -*-
r"""Test for PhaseTransition.py

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
# import shilofue.Foo as Foo  # import test module
from shilofue.PhaseTransition import *
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_parse_phase_input():
    '''
    Test parse phase input
    '''
    Inputs = {
        "rho0" : 3300.0,\
        "drho" : [520.0, 250.0, 670.0, 840.0, 560.0, 597.0, 1170.0],\
        "xc" : [0.55, 0.55, 0.05, 0.55, 0.4, 0.4, 0.4]
    }
    phase_opt = PHASE_OPT()
    phase_opt.import_options(Inputs)
    output = ParsePhaseInput(phase_opt)
    assert(output == "3300.0|3586.0|3723.5|3757.0|4219.0|4443.0|4681.8|5149.8")

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

