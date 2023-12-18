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
import filecmp  # for compare file contents
import numpy as np
from shilofue.Rheology import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')

# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_rheology')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_ConvertFCoh():
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    Operator = RHEOLOGY_OPR()
    # read profile
    # debug
    Operator.ReadProfile(da_file)
    rheology_aspect, _ = Operator.MantleRheology(rheology="HK03_wet_mod", debug=True)

    rheology = "HK03_wet_mod"
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)
    diffusion_creep_coh = ConvertFCoh(diffusion_creep)
    dislocation_creep_coh = ConvertFCoh(dislocation_creep)
    assert(abs(diffusion_creep_coh['A'] - 81787.10427350427) / 81787.10427350427 < 1e-6)
    assert(abs(dislocation_creep_coh['A'] - 2.1289006184840544) / 2.1289006184840544 < 1e-6)


def test_RheologyTableFormating():
    '''
    test function RheologyTableFormating
    assert:
        the latex table contents
    '''
    rheology = "HK03_f" 
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)
    header, data = RheologyTableFormating(diffusion_creep, "HK03 diffusion", format="HK03")
    TexTable = Utilities.TEX_TABLE("table-rheology",\
                                    header=header, data=data) # class initiation
    table_contents = TexTable(format="latex", fix_underscore_in_content=False)
    file_out = os.path.join(test_dir, "rheology_latex_table.tex")
    with open(file_out, 'w') as fout:
        fout.write(table_contents)
    assert(os.path.isfile(file_out))
    file_out_std = os.path.join(source_dir, "rheology_latex_table.tex")
    assert(filecmp.cmp(file_out, file_out_std))


def test_rheology_json():
    '''
    Test the implementation of the class RHEOLOGY_JSON
    '''
    # make a directory for the testing
    test_output_dir = os.path.join(test_dir, "test_rheology_json")
    if os.path.isdir(test_output_dir):
        rmtree(test_output_dir)
    os.mkdir(test_output_dir)
    # read json file
    json_file = os.path.join(source_dir, "test_rheology.json")
    RheologyJson = RHEOLOGY_JSON()
    RheologyJson.read_json(json_file)
    RheologyJson.check()
    # test utilities for PlotShearZoneRheologySummaryJson
    fig_path = os.path.join(test_output_dir, "test_plot_rate_stress.png")
    rheologies = RheologyJson.GetRheologyFeatures()
    rheologyOpt = rheologies[0]
    fig, ax = plt.subplots()
    grain_size = 3000 # u m
    PlotStrainRateStress(*rheologyOpt.to_RheologyInputs(), ax=ax)
    fig.savefig(fig_path)
    print("Save figure: ", fig_path)
    assert(os.path.isfile(fig_path))
    # test utilities for PlotViscosityTemperature
    fig_path = os.path.join(test_output_dir, "test_plot_visc_T.png")
    rheologies = RheologyJson.GetRheologyFeatures()
    rheologyOpt = rheologies[0]
    fig, ax = plt.subplots()
    grain_size = 3000 # u m
    strain_rate = 1e-13 # /s
    PlotViscosityTemperature(*rheologyOpt.to_RheologyInputs(), ax=ax, strain_rate=strain_rate)
    fig.savefig(fig_path)
    print("Save figure: ", fig_path)
    assert(os.path.isfile(fig_path))
    json_file = os.path.join(source_dir, "test_rheology_summary.json")
    
    PlotRheologySummaryJson(json_file)


def test_HK03_mod_whole_mantle_apsect():
    """
    test converting the HK03 modified rheology:
        1. to a whole-mantle rheology including the upper mantle and lower mantle rheologies
        2. to inputs in aspect
    """
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(da_file)
    rheology_aspect, _ = Operator.MantleRheology(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=20e3,\
    dVdiff=-5.5e-6, dVdisl=-1.2e-6, save_profile=1, debug=True)
    diffusion_creep = rheology_aspect['diffusion_creep']
    diffusion_A = diffusion_creep['A']
    diffusion_E = diffusion_creep['E']
    diffusion_V = diffusion_creep['V']
    assert(abs(diffusion_A - 1.2268065641025644e-16)/1.2268065641025644e-16 < 1e-6)
    assert(abs(diffusion_E - 285000.0) / 285000.0 < 1e-6)
    assert(abs(diffusion_V - 6.9e-6)/6.9e-6 < 1e-6)
    dislocation_creep = rheology_aspect['dislocation_creep']
    dislocation_A = dislocation_creep['A']
    dislocation_E = dislocation_creep['E']
    dislocation_V = dislocation_creep['V']
    assert(abs(dislocation_A - 5.019358499849437e-17)/5.019358499849437e-17 < 1e-6)
    assert(abs(dislocation_E - 480000.0)/480000.0 < 1e-6)
    assert(abs(dislocation_V - 1.008e-05)/1.008e-05 < 1e-6)
    diffusion_creep_lm = rheology_aspect['diffusion_lm']
    diffusion_lm_A = diffusion_creep_lm['A']
    diffusion_lm_V = diffusion_creep_lm['V']
    print("diffusion_creep_lm: ", diffusion_creep_lm) 
    assert(abs(diffusion_lm_A - 1.3427402992734117e-19)/1.3427402992734117e-19 < 1e-6)
    assert(abs(diffusion_lm_V - 3e-6)/3e-6 < 1e-6)

def test_HK03_mod_whole_mantle_apsect_TwoDSubdution():
    """
    test converting the HK03 modified rheology to a whole-mantle rheology in aspect
    and converts to a prm file
    Here I will use the MantleRheology where I tried to take care of the prefactor F.
    And I tried to fit the ones I have in v0 by modify the prefactors. 
    In this way I could document the form of flow law I used there (as I don't want to change that).
    Assertion:
        I could reproduce the parameters I have for the TwoDSubduction case.
    """
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(da_file)
    rheology_aspect, _ = Operator.MantleRheology(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=30e3,\
    dVdiff=-5.5e-6, dVdisl=2.12e-6, save_profile=1, dAdiff_ratio=0.33333247873, dAdisl_ratio=1.040297619,\
    jump_lower_mantle=15.0, debug=True)
    diffusion_creep = rheology_aspect['diffusion_creep']
    diffusion_A = diffusion_creep['A']
    diffusion_E = diffusion_creep['E']
    diffusion_V = diffusion_creep['V']
    assert(abs(diffusion_A - 4.089344729345424e-17)/4.089344729345424e-17 < 1e-6)
    assert(abs(diffusion_E - 285000.0) / 285000.0 < 1e-6)
    assert(abs(diffusion_V - 6.9e-6)/6.9e-6 < 1e-6)
    dislocation_creep = rheology_aspect['dislocation_creep']
    dislocation_A = dislocation_creep['A']
    dislocation_E = dislocation_creep['E']
    dislocation_V = dislocation_creep['V']
    assert(abs(dislocation_A - 5.2216266963007814e-17)/5.2216266963007814e-17 < 1e-6)
    assert(abs(dislocation_E - 490000.0)/490000.0 < 1e-6)
    assert(abs(dislocation_V - 1.34e-5)/1.34e-5 < 1e-6)
    diffusion_creep_lm = rheology_aspect['diffusion_lm']
    diffusion_lm_A = diffusion_creep_lm['A']
    diffusion_lm_V = diffusion_creep_lm['V']
    assert(abs(diffusion_lm_A - 6.250893063953478e-20)/6.250893063953478e-20 < 1e-6)
    assert(abs(diffusion_lm_V - 3e-6)/3e-6 < 1e-6)


def test_AB17_wet_whole_mantle_apsect_prm():
    """
    test converting the AB17 rheology to a whole-mantle rheology in aspect
    and converts to a prm file
    Here I will use the MantleRheology where I tried to take care of the prefactor F
    correctly
    """
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(da_file)
    rheology_aspect, _ = Operator.MantleRheology(rheology="AB17", save_profile=1)
    diffusion_creep = rheology_aspect['diffusion_creep']
    diffusion_A = diffusion_creep['A']
    diffusion_E = diffusion_creep['E']
    diffusion_V = diffusion_creep['V']
    # assert(abs(diffusion_A - 3.0000000000000002e-15)/3.0000000000000002e-15 < 1e-6)
    # assert(abs(diffusion_E - 335000.0) /335000.0  < 1e-6)
    # assert(abs(diffusion_V - 4e-6)/4e-6 < 1e-6)
    dislocation_creep = rheology_aspect['dislocation_creep']
    dislocation_A = dislocation_creep['A']
    dislocation_E = dislocation_creep['E']
    dislocation_V = dislocation_creep['V']
    # assert(abs(dislocation_A - 2.4007134284958627e-14)/2.4007134284958627e-14 < 1e-6)
    # assert(abs(dislocation_E - 480000.0)/480000.0 < 1e-6)
    # assert(abs(dislocation_V - 1.1e-5)/1.1e-5 < 1e-6)
    diffusion_creep_lm = rheology_aspect['diffusion_lm']
    diffusion_lm_A = diffusion_creep_lm['A']
    diffusion_lm_V = diffusion_creep_lm['V']
    # assert(abs(diffusion_lm_A - 1.0225822662706545e-16)/1.0225822662706545e-16 < 1e-6)
    # assert(abs(diffusion_lm_V - 3e-6)/3e-6 < 1e-6)
    pass


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

