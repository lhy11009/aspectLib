# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from copy import deepcopy
import shilofue.Cases as CasesP
import shilofue.ParsePrm as ParsePrm
from shilofue.Rheology import RHEOLOGY_OPR, RHEOLOGY_PRM, ConvertFromAspectInput,\
                        STRENGTH_PROFILE, GetRheology, Convert2AspectInput,\
                        AssignAspectViscoPlasticPhaseRheology, RefitRheology
from shilofue.TwoDSubduction0.Cases import re_write_geometry_while_assigning_plate_age, change_field_to_particle
import shilofue.Rheology_old_Dec_2023 as Rheology_old_Dec_2023

import json, re

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

year = 365 * 24 * 3600.0  # yr to s

twod_default_wb_file = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "220716", "case.wb")

class CASE_OPT(CasesP.CASE_OPT):
    '''
    Define a class to work with CASE
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        see document (run with -h option) for detail
        '''
        CasesP.CASE_OPT.__init__(self)
        self.start = self.number_of_keys()
        self.add_key("Width of the box", float, ['geometry setup', 'box width'], 2000e3, nick='box_width')
        self.add_key("Length of the box", float, ['geometry setup', 'box length'], 4000e3, nick='box_length')
        self.add_key("Height of the box", float, ['geometry setup', 'box height'], 1000e3, nick='box_height')
        self.add_key("Width of the subducting plate", float, ['plate setup', 'sp width'], 1000e3, nick='sp_width')
        self.add_key("Length of the subducting plate", float, ['plate setup', 'sp length'], 2200e3, nick='sp_length')
        self.add_key("Depth of the subducting plate (in schellart 07 model", float,\
                    ['plate setup', 'sp depth'], 100e3, nick='sp_depth')
        self.add_key("Length of the trailing tail of the plate", float, ['plate setup', 'trailing length'], 275e3, nick='trailing_length')
        self.add_key("Reset composition for the trailing tail of the plate", int, ['plate setup', 'reset trailing morb'], 0, nick='reset_trailing_morb')
        self.add_key("Reference viscosity of the mantle", float, ['rheology', 'reference viscosity'], 1e20, nick='ref_visc')
        self.add_key("Relative viscosity of the plate", float, ['rheology', 'plate relative viscosity'], 200.0, nick='relative_visc_plate')
        self.add_key("Yielding friction angle", float, ['rheology', 'friction angle'], 5.71, nick='friction_angle')
        self.add_key("Relative viscosity of the lower mantle", float, ['rheology', 'lower mantle relative viscosity'], 100.0, nick='relative_visc_lower_mantle')
        self.add_key("Yielding cohesion", float, ['rheology', 'cohesion'], 48e6, nick='cohesion')
        self.add_key("Depth of the subducting plate to refilen (in schellart 07 model)", float,\
                    ['plate setup', 'sp depth refining'], 200e3, nick='sp_depth_refining')
        self.add_key("Mantle reference density", float, ['reference density'], 3300.0, nick='reference_density')
        self.add_key("Excess of density of the subducting plate (in schellart 07 model)", float,\
                    ['plate setup', 'sp relative density'], 80.0, nick='sp_relative_density')
        self.add_key("Global refinement", int, ['refinement', 'global refinement'], 4, nick='global_refinement')
        self.add_key("Adaptive refinement", int, ['refinement', 'adaptive refinement'], 2, nick='adaptive_refinement')
        self.add_key("mantle rheology", str, ['mantle rheology', 'scheme'], "HK03_wet_mod_twod", nick='mantle_rheology_scheme')
        self.add_key("Thickness of the shear zone / crust", float, ["shear zone", 'thickness'], 50e3, nick='Dsz')
        self.add_key("Thickness of the depleted lithosphere", float, ['plate setup', 'dl thickness'], 50e3, nick='Ddl')
        self.add_key("Apply the mantle reference density for all the compositions", int, ['apply reference density'],\
        0, nick='apply_reference_density')
        self.add_key("Length of the initial slab", float, ['slab setup', 'length'], 167e3, nick='slab_length')
        self.add_key("Dipping angle of the initial slab", float, ['slab setup', 'dip'], 15.5, nick='dip_angle')
        self.add_key("Age of the subducting plate at trench", float, ['plate setup', 'sp age'], 80e6, nick='sp_age_trench')
        self.add_key("Age of the overiding plate", float, ['plate setup', 'ov age'], 40e6, nick='ov_age')
        self.add_key("Method to generate the geometry", str, ['setup method'],'manual', nick='setup_method')
        self.add_key("Spreading velocity of the sp plate, use with the method \"2d_consistent\"", float, ['plate setup', 'sp rate'], 0.05, nick='sp_rate')
        self.add_key("Length of the Box before adjusting for the age of the trench.\
This is used with the option \"adjust box width\" for configuring plate age at the trench.\
This value is the width of box for a default age (i.e. 80Myr), while the width of box for a\
different age will be adjusted.",\
          float, ["geometry setup", "box length before adjusting"], 6.783e6, nick='box_length_pre_adjust')
        self.add_key("Reset viscosity for the trailing tail of the overiding plate", int, ['rheology', 'reset trailing ov viscosity'], 0, nick='reset_trailing_ov_viscosity')
        self.add_key("viscous flow law for mantle rheology", str, ['mantle rheology', 'flow law'], "diffusion", nick='mantle_rheology_flow_law')
        self.add_key("use WB new ridge implementation", int, ['world builder', 'use new ridge implementation'], 0, nick='wb_new_ridge')
        self.add_key("Assign a side plate", int, ['plate setup', 'assign side plate'], 0, nick='assign_side_plate')
        self.add_key("branch", str, ['branch'], "", nick='branch')
        self.add_key("Age of the transit overiding plate", float, ['plate setup', 'ov transit age'], -1.0, nick='ov_trans_age')
        self.add_key("Length of the transit overiding plate", float,\
            ['plate setup', 'ov transit length'], 300e3, nick='ov_trans_length')
        self.add_key("Geometry", str, ["geometry"], "box", nick='geometry')
        self.add_key("position of the sp ridge", float, ['plate setup', 'sp ridge x'], 0.0, nick='sp_ridge_x')
        self.add_key("distance of the ov ridge from the side wall", float, ['plate setup', 'ov side dist'], 0.0, nick='ov_side_dist')
        self.add_key("prescribe mantle temperature before the subducting plate starts",\
        int, ['plate setup', 'prescribe mantle sp start'], 0, nick='prescribe_mantle_sp')
        self.add_key("prescribe mantle temperature after the overiding plate ends",\
        int, ['plate setup', 'prescribe mantle ov end'], 0, nick='prescribe_mantle_ov')
        self.add_key("Minimum mantle rheology for subducting initiation",\
        float, ['mantle rheology', 'minimum viscosity for initial stage'], -1.0, nick='mantle_minimum_init')
        self.add_key("prescribe viscosity for composition",\
        int, ['plate setup', 'reset composition viscosity'], 0, nick='reset_composition_viscosity')
        self.add_key("prescribe viscosity for composition width",\
        float, ['plate setup', 'reset composition viscosity width'], 2000e3, nick='reset_composition_viscosity_width')
        self.add_key("adjust the length of the box by adding the trailing length to the box",\
        int, ['geometry setup', 'adjust box trailing length'], 0, nick='adjust_box_trailing_length')
        self.add_key("method to set up the slices",\
        str, ['geometry setup', 'repitition slice method'], 'floor', nick='repitition_slice_method')
        self.add_key("Viscosity in the slab core", float,\
            ['shear zone', "slab core viscosity"], -1.0, nick='slab_core_viscosity')
        self.add_key("Minimum viscosity", float,\
            ["minimum viscosity"], 1e19, nick='global_minimum_viscosity')
        self.add_key("Make 2d consistent plate, only works with the 2d_consistent setup", int,\
            ["make 2d consistent plate"], 0, nick='make_2d_consistent_plate')
        self.add_key("If the side of the box is coarsened in assigning the minimum refinement function", int,\
            ["refinement", "coarsen side"], 0, nick='coarsen_side')
        self.add_key("The side of the box is coarsened, except for an interval attached to the plate side", float,\
            ["refinement", "coarsen side interval"], -1.0, nick='coarsen_side_interval')
        self.add_key("automatically fix boundary temperature",\
        int, ['geometry setup', 'fix boudnary temperature auto'], 0, nick='fix_boudnary_temperature_auto')
        self.add_key("The side of the box is coarsened with this level", int,\
            ["refinement", "coarsen side level"], -1, nick='coarsen_side_level')
        self.add_key("The difference between the highest level used in the MR function and the total refinement level",\
         int, ["refinement", "coarsen minimum refinement level"], 1, nick='coarsen_minimum_refinement_level')
        self.add_key("Include peierls creep", int, ['include peierls creep'], 0, nick='if_peierls')
        self.add_key("Fix the activation volume of the Peierls creep", str,\
         ['peierls creep', "fix peierls V as"], '', nick="fix_peierls_V_as")
        self.add_key("Length of the trailing tail of the plate, another implementation", float, ['plate setup', 'trailing length 1'],\
                     0.0, nick='trailing_length_1')
        self.add_key("Value of Coh to use in the rheology", float,\
            ['mantle rheology', "Coh"], 500.0, nick='mantle_coh')
        self.add_key("Value of lower/upper mantle ratio to use in the rheology", float,\
            ['mantle rheology', "jump lower mantle"], 100.0, nick='jump_lower_mantle')
        self.add_key("Adjust detail of mantle rheology", int,\
            ['mantle rheology', "adjust detail"], 0, nick='adjust_mantle_rheology_detail')
        self.add_key("Include upper plate composition for the overriding plate",\
        int, ['plate setup', 'include ov upper plate'], 0, nick='include_ov_upper_plate')
        # todo_strength
        self.add_key("Slab strengh", float, ['plate setup', "strength"], 500e6, nick="slab_strength")

    
    def check(self):
        _type = self.values[9] 
        reset_trailing_morb = self.values[self.start+7]
        assert(reset_trailing_morb in [0, 1, 2, 3])
        friction_angle = self.values[self.start+10] # range of friction angle, in degree
        assert(friction_angle >= 0.0 and friction_angle <= 90.0)
        dip_angle = self.values[self.start+23] # initial dipping angle of the slab
        assert(dip_angle > 0 and dip_angle <= 90.0) # an angle between 0 and 90
        setup_method = self.values[self.start+26] # method of seting up slabs
        assert(setup_method in ['manual', '2d_consistent'])
        repitition_slice_method = self.values[self.start+45]
        assert(repitition_slice_method in ['floor', 'nearest'])
        make_2d_consistent_plate = self.values[self.start+48]
        coarsen_side = self.values[self.start+49]
        if make_2d_consistent_plate > 0:
            assert(_type == "2d_consistent" and setup_method == '2d_consistent')
        if coarsen_side > 0:
            assert(_type == "2d_consistent")

    def reset_refinement(self, refinement_level):
        '''
        Note:
            reload function from base class
        '''
        self.values[self.start+16] = refinement_level

    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        if_wb = self.values[8]
        _type = self.values[9] 
        setup_method = self.values[self.start+26] # method of seting up slabs
        box_width = self.values[self.start]
        box_length = self.values[self.start+1]
        if setup_method == '2d_consistent':
            # adjust box width with the age and plate spreading rate
            box_length = re_write_geometry_while_assigning_plate_age(
            *self.to_re_write_geometry_pa()
            ) # adjust box width
        box_height = self.values[self.start+2]
        sp_width = self.values[self.start+3]
        trailing_length = self.values[self.start+6]
        reset_trailing_morb = self.values[self.start+7]
        ref_visc = self.values[self.start+8]
        relative_visc_plate = self.values[self.start+9]
        friction_angle = self.values[self.start+10]
        relative_visc_lower_mantle = self.values[self.start+11]
        cohesion = self.values[self.start+12]
        sp_depth_refining = self.values[self.start+13]
        reference_density = self.values[self.start+14]
        sp_relative_density = self.values[self.start+15]
        global_refinement = self.values[self.start+16]
        adaptive_refinement = self.values[self.start+17]
        mantle_rheology_scheme = self.values[self.start+18]
        Dsz = self.values[self.start+19]
        Ddl = self.values[self.start+20]
        apply_reference_density = self.values[self.start+21]
        reset_trailing_ov_viscosity = self.values[self.start+29]
        mantle_rheology_flow_law = self.values[self.start+30]
        stokes_solver_type = self.values[18]
        case_o_dir = self.values[16]
        branch = self.values[self.start+33]
        geometry = self.values[self.start+36]
        sp_ridge_x = self.values[self.start+37]
        ov_side_dist = self.values[self.start+38]
        prescribe_mantle_sp = self.values[self.start+39]
        prescribe_mantle_ov = self.values[self.start+40]
        mantle_minimum_init = self.values[self.start+41]
        visual_software = self.values[24] 
        comp_method = self.values[25] 
        reset_composition_viscosity = self.values[self.start+42]
        reset_composition_viscosity_width = self.values[self.start+43]
        slab_core_viscosity = self.values[self.start+46]
        global_minimum_viscosity = self.values[self.start+47]
        if visual_software == 'paraview':
            # output the step 1 if the fast_first_step is processed
            self.output_step_one_with_fast_first_step()
        repitition_slice_method = self.values[self.start+45]
        coarsen_side = self.values[self.start+49]
        coarsen_side_interval = self.values[self.start+50]
        fix_boudnary_temperature_auto = self.values[self.start+51]
        coarsen_side_level = self.values[self.start+52]
        coarsen_minimum_refinement_level = self.values[self.start+53]
        use_new_rheology_module = self.values[28]
        if_peierls = self.values[self.start+54]
        fix_peierls_V_as = self.values[self.start+55]
        trailing_length_1 = self.values[self.start+56]
        sp_rate = self.values[self.start+27] # method of seting up slabs
        ov_age = self.values[self.start+25]
        detail_mantle_coh = self.values[self.start+57]
        detail_jump_lower_mantle = self.values[self.start+58]
        adjust_mantle_rheology_detail = self.values[self.start+59]
        include_ov_upper_plate = self.values[self.start+60]
        slab_strength = self.values[self.start+61]
        return _type, if_wb, geometry, box_width, box_length, box_height,\
            sp_width, trailing_length, reset_trailing_morb, ref_visc,\
            relative_visc_plate, friction_angle, relative_visc_lower_mantle, cohesion,\
            sp_depth_refining, reference_density, sp_relative_density, global_refinement,\
            adaptive_refinement, mantle_rheology_scheme, Dsz, apply_reference_density, Ddl,\
            reset_trailing_ov_viscosity, mantle_rheology_flow_law, stokes_solver_type, case_o_dir,\
            branch, sp_ridge_x, ov_side_dist, prescribe_mantle_sp, prescribe_mantle_ov, mantle_minimum_init,\
            comp_method, reset_composition_viscosity, reset_composition_viscosity_width, repitition_slice_method,\
            slab_core_viscosity, global_minimum_viscosity, coarsen_side, coarsen_side_interval, fix_boudnary_temperature_auto,\
            coarsen_side_level, coarsen_minimum_refinement_level, use_new_rheology_module, if_peierls, fix_peierls_V_as,\
            trailing_length_1, sp_rate, ov_age, detail_mantle_coh, detail_jump_lower_mantle, adjust_mantle_rheology_detail,\
            include_ov_upper_plate, slab_strength
        
    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[8]
        _type = self.values[9] 
        sp_width = self.values[self.start+3]
        sp_length = self.values[self.start+4]
        trailing_length = self.values[self.start+6]
        Dsz = self.values[self.start+19]
        Ddl = self.values[self.start+20]
        slab_length = self.values[self.start+22]
        dip_angle = self.values[self.start+23]
        sp_age_trench = self.values[self.start+24]
        ov_age = self.values[self.start+25]
        setup_method = self.values[self.start+26] # method of seting up slabs
        sp_rate = self.values[self.start+27] # method of seting up slabs
        wb_new_ridge = self.values[self.start+31]
        assign_side_plate = self.values[self.start+32]
        ov_trans_age =  self.values[self.start+34]
        ov_trans_length =  self.values[self.start+35]
        if ov_trans_age < 0.0:
            if_ov_trans = False
        else:
            if_ov_trans = True
        geometry = self.values[self.start+36]
        sp_ridge_x = self.values[self.start+37]
        ov_side_dist = self.values[self.start+38]
        box_length = self.values[self.start+1]
        if setup_method == '2d_consistent':
            # adjust box width with the age and plate spreading rate
            box_length = re_write_geometry_while_assigning_plate_age(
            *self.to_re_write_geometry_pa()
            ) # adjust box width
        make_2d_consistent_plate = self.values[self.start+48]
        return _type, if_wb, geometry, sp_width, sp_length, trailing_length,\
        Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age,\
        setup_method, sp_rate, wb_new_ridge, assign_side_plate,\
        if_ov_trans, ov_trans_age, ov_trans_length, sp_ridge_x,\
        ov_side_dist, box_length, make_2d_consistent_plate
    
    def to_re_write_geometry_pa(self):
        '''
        Interface to re_write_geometry_while_assigning_plate_age
        '''
        box_length_pre_adjust = self.values[self.start+28]
        default_sp_age_trench = self.defaults[self.start+24]
        sp_age_trench = self.values[self.start+24]
        sp_rate = self.values[self.start+27] # method of seting up slabs
        trailing_length = self.values[self.start+6]
        # todo_adjust
        adjust_box_trailing_length = self.values[self.start+44]
        if adjust_box_trailing_length:
            adjust_trailing_length = trailing_length
        else:
            adjust_trailing_length = 0.0
        # the 0.0 is appended for the sp_trailing_length
        return box_length_pre_adjust, default_sp_age_trench, sp_age_trench, sp_rate, adjust_trailing_length, 0.0


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, _type, if_wb, geometry, box_width, box_length, box_height,\
    sp_width, trailing_length, reset_trailing_morb, ref_visc, relative_visc_plate, friction_angle,\
    relative_visc_lower_mantle, cohesion, sp_depth_refining, reference_density, sp_relative_density, \
    global_refinement, adaptive_refinement, mantle_rheology_scheme, Dsz, apply_reference_density, Ddl,\
    reset_trailing_ov_viscosity, mantle_rheology_flow_law, stokes_solver_type, case_o_dir, branch,\
    sp_ridge_x, ov_side_dist, prescribe_mantle_sp, prescribe_mantle_ov, mantle_minimum_init, comp_method,\
    reset_composition_viscosity, reset_composition_viscosity_width, repitition_slice_method, slab_core_viscosity,\
    global_minimum_viscosity, coarsen_side, coarsen_side_interval, fix_boudnary_temperature_auto, coarsen_side_level,\
    coarsen_minimum_refinement_level, use_new_rheology_module, if_peierls, fix_peierls_V_as, trailing_length_1, sp_rate,\
    ov_age, detail_mantle_coh, detail_jump_lower_mantle, adjust_mantle_rheology_detail, include_ov_upper_plate, slab_strength):
        '''
        Configure prm file
        '''
        self.configure_case_output_dir(case_o_dir)
        o_dict = self.idict.copy()
        if branch != "":
            if branch == "master":
                branch_str = ""
            else:
                branch_str = "_%s" % branch
            o_dict["Additional shared libraries"] =  "$ASPECT_SOURCE_DIR/build%s/visco_plastic_TwoD/libvisco_plastic_TwoD.so" % branch_str
            if "Prescribe internal mantle adiabat temperatures" in o_dict:
                # append the library to prescribe temperature if needed
                o_dict["Additional shared libraries"] += ", "
                o_dict["Additional shared libraries"] +=  "$ASPECT_SOURCE_DIR/build%s/prescribe_field_T_adiabat/libprescribe_field_T_adiabat.so" % branch_str
            # another implementation, set trailing_length = 0.0 to turn off the older implementation
            # here we set this option to false and remove the corresponding section
            if trailing_length < 1e-6 and trailing_length_1 > 1e-6:
                o_dict["Additional shared libraries"] += ", "
                o_dict["Additional shared libraries"] +=  "$ASPECT_SOURCE_DIR/build%s/prescribe_field/libprescribed_temperature.so" % branch_str
        # geometry options
        # Box size: assigned 
        # repitition: figure this out by deviding the dimensions with a unit value of repitition_slice
        # The repitition_slice is defined by a max_repitition_slice and the minimum value compared with the dimensions
        # boudnary temperature: figure this out from the depth average profile
        max_repitition_slice = 1000e3 # a value for the maximum value
        repitition_slice = np.min(np.array([box_length, box_width, box_height, max_repitition_slice]))  # take the min as the repitition_slice
        Ro = 6371e3
        box_length_deg = box_length / Ro * 180.0 / np.pi
        box_width_deg = box_width / Ro * 180.0 / np.pi
        # invoke different options for different geometries
        if geometry == "box":
            o_dict['Geometry model']['Box']['X extent'] =  str(box_length)
            o_dict['Geometry model']['Box']['Y extent'] =  str(box_width)
            o_dict['Geometry model']['Box']['Z extent'] =  str(box_height)
            if repitition_slice_method == "floor": 
                o_dict['Geometry model']['Box']['X repetitions'] = str(int(box_length//repitition_slice))
                o_dict['Geometry model']['Box']['Y repetitions'] = str(int(box_width//repitition_slice))
                o_dict['Geometry model']['Box']['Z repetitions'] = str(int(box_height//repitition_slice))
            elif repitition_slice_method == "nearest": 
                o_dict['Geometry model']['Box']['X repetitions'] = \
                    str(int(np.ceil(int((box_length/repitition_slice) * 2.0) / 2.0)))
                o_dict['Geometry model']['Box']['Y repetitions'] = \
                    str(int(np.ceil(int((box_width/repitition_slice) * 2.0) / 2.0)))
                o_dict['Geometry model']['Box']['Z repetitions'] = \
                    str(int(np.ceil(int((box_height/repitition_slice) * 2.0) / 2.0)))
        elif geometry == "chunk":
            o_dict['Geometry model']['Model name'] = 'chunk'
            o_dict['Geometry model'].pop('Box')
            o_dict['Geometry model']['Chunk'] = {}
            o_dict['Geometry model']['Chunk']['Chunk inner radius'] =  "%.4e" % (Ro-box_height)
            o_dict['Geometry model']['Chunk']['Chunk outer radius'] = "6.371e6" 
            o_dict['Geometry model']['Chunk']["Chunk minimum longitude"] = "0"
            o_dict['Geometry model']['Chunk']['Chunk maximum longitude'] =  str(box_length_deg)
            o_dict['Geometry model']['Chunk']["Chunk minimum latitude"] = "0"
            o_dict['Geometry model']['Chunk']['Chunk maximum latitude'] =  str(box_width_deg)
            if repitition_slice_method == "floor": 
                o_dict['Geometry model']['Chunk']['Longitude repetitions'] = str(int(box_length//repitition_slice))
                o_dict['Geometry model']['Chunk']['Latitude repetitions'] = str(int(box_width//repitition_slice))
                o_dict['Geometry model']['Chunk']['Radius repetitions'] = str(int(box_height//repitition_slice))
            elif repitition_slice_method == "nearest": 
                o_dict['Geometry model']['Chunk']['Longitude repetitions'] = \
                    str(int(np.ceil(int((box_length/repitition_slice) * 2.0) / 2.0)))
                o_dict['Geometry model']['Chunk']['Latitude repetitions'] = \
                    str(int(np.ceil(int((box_width/repitition_slice) * 2.0) / 2.0)))
                o_dict['Geometry model']['Chunk']['Radius repetitions'] = \
                    str(int(np.ceil(int((box_height/repitition_slice) * 2.0) / 2.0)))

        # boundary temperature
        if fix_boudnary_temperature_auto:
            # boudnary temperature: figure this out from the depth average profile
            assert(self.da_Tad_func is not None)
            try:
                Tad_bot = self.da_Tad_func(box_height) # bottom adiabatic temperature
            except ValueError:
                # in case this is above the given range of depth in the depth_average
                # file, apply a slight variation and try again
                Tad_bot = self.da_Tad_func(box_height - 50e3)
            o_dict['Boundary temperature model']['Box']['Bottom temperature'] = "%.4e" % Tad_bot
        # modify for the chunk geometry
        if geometry == "chunk":
            o_dict['Boundary temperature model']["List of model names"] = "spherical constant"
            o_dict['Boundary temperature model']["Spherical constant"] = o_dict['Boundary temperature model'].pop("Box")
            o_dict['Boundary temperature model']["Spherical constant"]["Inner temperature"] = o_dict['Boundary temperature model']["Spherical constant"].pop("Bottom temperature")
            o_dict['Boundary temperature model']["Spherical constant"]["Outer temperature"] = o_dict['Boundary temperature model']["Spherical constant"].pop("Top temperature")

        # refinement
        o_dict["Mesh refinement"]["Initial global refinement"] = str(global_refinement)
        o_dict["Mesh refinement"]["Minimum refinement level"] = str(global_refinement)
        o_dict["Mesh refinement"]["Initial adaptive refinement"] = str(adaptive_refinement)
        # Minimum refinement function
        # the largest value used in the minimum refinement function is determined by
        # max_refinement - coarsen_minimum_refinement_level
        max_refinement = global_refinement + adaptive_refinement + 1
        if _type == "s07":
            if (abs(sp_depth_refining - 200e3)/200e3 > 1e-6):
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3" % (box_height, sp_depth_refining)
            else:
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=200e3, Dp=100e3" % (box_height)
        elif _type == "s07T":
            o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3, Rd=%d, Rum=%d" % (box_height, sp_depth_refining,\
                max_refinement - coarsen_minimum_refinement_level, max_refinement - coarsen_minimum_refinement_level - 1)
        elif _type == "2d_consistent":
            if geometry == "box":
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                    "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3, Rd=%d, Rum=%d" % (box_height, sp_depth_refining,\
                    max_refinement - coarsen_minimum_refinement_level, max_refinement - coarsen_minimum_refinement_level -1)
            if geometry == "chunk":
                o_dict["Mesh refinement"]["Minimum refinement function"]["Coordinate system"] = "spherical"
                o_dict["Mesh refinement"]["Minimum refinement function"]["Variable names"] = "r,phi,theta"
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                    "Ro=%.4e, UM=670e3, DD=%.4e, Rd=%d, Rum=%d" % (Ro, sp_depth_refining,\
                    max_refinement - coarsen_minimum_refinement_level, max_refinement - coarsen_minimum_refinement_level -1)
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function expression"] = "((Ro-r<UM)? ((Ro-r<DD)? Rd: Rum): 0.0)"
            if coarsen_side:
                if geometry == "chunk":
                    raise NotImplementedError()
                # append an additional variable to coarsen the side with no slab
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                    "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3, Wside=%.4e, Rd=%d, Rum=%d" % (box_height, sp_depth_refining,\
                    sp_width + coarsen_side_interval,\
                    max_refinement - coarsen_minimum_refinement_level, max_refinement - coarsen_minimum_refinement_level - 1)
                if coarsen_side_level ==-1:
                    # in this case, coarsen all
                    o_dict["Mesh refinement"]["Minimum refinement function"]["Function expression"]=\
                            " (y < Wside)?\\\n\
                                    ((Do-z<UM)?\\\n\
                                      ((Do-z<DD)?\\\n\
				                        ((Do-z<Dp+50e3)? Rd: Rum)\\\n\
				                        :Rum)\\\n\
                                      :0)\\\n\
                                    :0"
                else:
                    # in this case, coarsen by a certain level
                    o_dict["Mesh refinement"]["Minimum refinement function"]["Function expression"]=\
                            " (y < Wside)?\\\n\
                                    ((Do-z<UM)?\\\n\
                                      ((Do-z<DD)?\\\n\
				                        ((Do-z<Dp+50e3)? Rd: Rum)\\\n\
				                        :Rum)\\\n\
                                      :0)\\\n\
				                    :((Do-z<UM)?\\\n\
                                      ((Do-z<DD)?\\\n\
				                        ((Do-z<Dp+50e3)? Rd-%d: Rum-%d)\\\n\
				                        :Rum-%d)\\\n\
                                      :0)" % (coarsen_side_level, coarsen_side_level, coarsen_side_level)
                    pass
        # composition method
        if comp_method == "field":
            pass
        elif comp_method == "particle":
            o_dict = change_field_to_particle(o_dict)
        
        # boundary temperature model
        # o_dict['Boundary temperature model'] =
        
        # Adiabatic surface temperature
        # o_dict["Adiabatic surface temperature"] = str(potential_T)
        # o_dict["Boundary temperature model"]["Box"]["Top temperature"] = str(potential_T)

        # materical model
        # a. modify material model
        # b. then merge it into the material model of the prm file.
        outputs = {}
        material_model_subsection = "Visco Plastic TwoD"
        # densities
        sp_density = reference_density + sp_relative_density
        if _type == 's07':
           if (abs(reference_density-3300.0)/3300.0 > 1e-6 or\
            abs(sp_density-3380.0)/3380.0 > 1e-6):
               o_dict['Material model'][material_model_subsection]['Densities'] =\
                   'background: %.4e, sp_upper: %.4e, sp_lower: %.4e'\
                   % (reference_density, sp_density, sp_density)
        if apply_reference_density:
            # apply reference density, rewrting all the previous settings
            o_dict['Material model'][material_model_subsection]['Densities'] = str(reference_density)
        # rheology
        # 1. mantle rheology
        if _type == 's07_newton':
            # for the type of model modified from schellart 07 paper and includes the newtonian
            # rheology. A profile is read from a depth_average file to get a reference temperature
            # and pressure.
            # Here I inherit the older rheology from my 2-d models:
            #   HK03: the original rheology from the Hirth and Kolstedt model
            #   HK03_wet_mod: the wet modified rheology I used in my 2-d models
            #   else: something else
            # use_new_rheology_module = 0 uses the old module. It's wrong.
            #   this option is left open to maintain backward compatibility
            da_file = os.path.join(ASPECT_LAB_DIR, 'files', 'ThDSubduction', "depth_average.txt")
            assert(os.path.isfile(da_file))
            if use_new_rheology_module == 1:
                Operator = RHEOLOGY_OPR()
            else:
                Operator = Rheology_old_Dec_2023.RHEOLOGY_OPR()
            # read profile
            Operator.ReadProfile(da_file)
            if mantle_rheology_scheme == "HK03":
                rheology, _ = Operator.MantleRheology(rheology="HK03",\
                            save_profile=1, save_json=1, use_effective_strain_rate=False)
            elif mantle_rheology_scheme == "HK03_wet_mod_twod":
                rheology, _ = Operator.MantleRheology(rheology="HK03_wet_mod",\
                            save_profile=1, save_json=1, use_effective_strain_rate=True,\
                            dEdiff=-40e3, dEdisl=30e3, dVdiff=-5.5e-6, dVdisl=2.12e-6,\
                            dAdiff_ratio=0.33333333333, dAdisl_ratio=1.040297619,\
                            jump_lower_mantle=15.0)
            elif mantle_rheology_scheme == "HK03_WarrenHansen23":
                # read in the original rheology
                rheology_name = "WarrenHansen23"
                rheology_prm_dict = RHEOLOGY_PRM()
                diffusion_creep_ori = getattr(rheology_prm_dict, rheology_name + "_diff")
                dislocation_creep_ori = getattr(rheology_prm_dict, rheology_name + "_disl")
                rheology_dict = {'diffusion': diffusion_creep_ori, 'dislocation': dislocation_creep_ori}
                # prescribe the correction
                diff_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': -2.1e-6}
                disl_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': 3e-6}
                # prescribe the reference state
                ref_state = {}
                if adjust_mantle_rheology_detail:
                    mantle_coh = detail_mantle_coh
                    jump_lower_mantle = detail_jump_lower_mantle
                else:
                    mantle_coh = 500.0
                    jump_lower_mantle = 100.0
                ref_state["Coh"] = mantle_coh # H / 10^6 Si
                ref_state["stress"] = 50.0 # MPa
                ref_state["P"] = 100.0e6 # Pa
                ref_state["T"] = 1250.0 + 273.15 # K
                ref_state["d"] = 15.0 # mu m
                # refit rheology
                rheology_dict_refit = RefitRheology(rheology_dict, diff_correction, disl_correction, ref_state)
                # derive mantle rheology
                rheology, _ = Operator.MantleRheology(assign_rheology=True, diffusion_creep=rheology_dict_refit['diffusion'],\
                                                            dislocation_creep=rheology_dict_refit['dislocation'], save_profile=1,\
                                                            use_effective_strain_rate=True, save_json=1, Coh=mantle_coh,\
                                                            jump_lower_mantle=jump_lower_mantle)
            else:
                rheology, _ = Operator.MantleRheology(rheology=mantle_rheology_scheme,\
                            save_profile=1, save_json=1, use_effective_strain_rate=True)
            s07T_assign_mantle_rheology(o_dict, rheology)    
            self.output_files.append(Operator.output_json)
            self.output_files.append(Operator.output_json_aspect)
            self.output_imgs.append(Operator.output_profile) # append plot of initial conition to figures
            # twick the type of flow law
            if mantle_rheology_flow_law == "composite":
                o_dict['Material model'][material_model_subsection]["Viscous flow law"] = "composite"
                o_dict = CasesP.SetNewtonSolver(o_dict)
        if _type == '2d_consistent':
            # for the 2d_consistent type
            # The adjustment of rheology is done after I implement the "HK03_WarrenHansen23" type
            da_file = os.path.join(ASPECT_LAB_DIR, 'files', 'ThDSubduction', "depth_average.txt")
            assert(os.path.isfile(da_file))
            Operator = RHEOLOGY_OPR()
            # read profile
            Operator.ReadProfile(da_file)
            if mantle_rheology_scheme == "HK03_WarrenHansen23":
                # read in the original rheology
                mantle_coh = 500 # be default, use a value of 500
                rheology_name = "WarrenHansen23"
                rheology_prm_dict = RHEOLOGY_PRM()
                diffusion_creep_ori = getattr(rheology_prm_dict, rheology_name + "_diff")
                dislocation_creep_ori = getattr(rheology_prm_dict, rheology_name + "_disl")
                rheology_dict = {'diffusion': diffusion_creep_ori, 'dislocation': dislocation_creep_ori}
                # prescribe the correction
                diff_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': -2.1e-6}
                disl_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': 3e-6}
                # prescribe the reference state
                if adjust_mantle_rheology_detail:
                    mantle_coh = detail_mantle_coh
                    jump_lower_mantle = detail_jump_lower_mantle
                else:
                    mantle_coh = 500.0
                    jump_lower_mantle = 100.0
                ref_state = {}
                ref_state["Coh"] = mantle_coh # H / 10^6 Si
                ref_state["stress"] = 50.0 # MPa
                ref_state["P"] = 100.0e6 # Pa
                ref_state["T"] = 1250.0 + 273.15 # K
                ref_state["d"] = 15.0 # mu m
                # refit rheology
                rheology_dict_refit = RefitRheology(rheology_dict, diff_correction, disl_correction, ref_state)
                # derive mantle rheology
                rheology, _ = Operator.MantleRheology(assign_rheology=True, diffusion_creep=rheology_dict_refit['diffusion'],\
                                                            dislocation_creep=rheology_dict_refit['dislocation'], save_profile=1,\
                                                            use_effective_strain_rate=True, save_json=1, Coh=mantle_coh,\
                                                            jump_lower_mantle=jump_lower_mantle)
            if mantle_rheology_scheme in ["HK03_WarrenHansen23"]: 
                consistent_2d_assign_mantle_rheology(o_dict, rheology, include_ov_upper_plate=include_ov_upper_plate)    
                self.output_files.append(Operator.output_json)
                self.output_files.append(Operator.output_json_aspect)
                self.output_imgs.append(Operator.output_profile) # append plot of initial conition to figures
            # Change the slab strength by the maximum yield stress, the default value is 500 Mpa
            # This implementation has a historical reason. I mistook 500e6 as the default value in aspect, but the 
            # actual default value is infinite (1e12). This is later adapted to using a large value.
            if abs(slab_strength - 500e6) / 500e6 > 1e-6 and slab_strength < 1e10:
                o_dict['Material model']['Visco Plastic TwoD']["Maximum yield stress"] = "%.4e" % slab_strength
            

        # 2. Yielding criteria
        if _type == 's07':
            prefactor_ref = 1.0 / 2.0 / ref_visc  # prefactors for diffusion creep
            prefactor_plate = 1.0 / 2.0 / (ref_visc * relative_visc_plate)
            prefactor_lower = 1.0 / 2.0 / (ref_visc * relative_visc_lower_mantle)
            o_dict['Material model'][material_model_subsection]['Prefactors for diffusion creep'] =\
                'background: %.4e|%.4e, sp_upper: %.4e, sp_lower: %.4e'\
                % (prefactor_ref, prefactor_lower, prefactor_plate, prefactor_plate)
            o_dict['Material model'][material_model_subsection]['Angles of internal friction'] =\
                "background:0.0, sp_upper: %.4e, sp_lower: 0.0" % friction_angle
            o_dict['Material model'][material_model_subsection]['Cohesions'] = "background:1e31, sp_upper: %.4e, sp_lower:1e31" % cohesion
        # 3. set the minium rheology
        if mantle_minimum_init > 0.0:
            o_dict['Material model'][material_model_subsection]['Minimum viscosity'] = "%.4e" % mantle_minimum_init
        if _type == "2d_consistent" and slab_core_viscosity > 0.0:
            # assign a strong core inside the slab
            o_dict['Material model'][material_model_subsection]['Minimum viscosity'] =\
                "background: %.4e, sp_upper: %.4e, sp_lower: %.4e, plate_edge: %.4e" %\
                    (global_minimum_viscosity, global_minimum_viscosity, slab_core_viscosity, global_minimum_viscosity)
            if include_ov_upper_plate:
                o_dict['Material model'][material_model_subsection]['Minimum viscosity'] =\
                    "background: %.4e, sp_upper: %.4e, sp_lower: %.4e, plate_edge: %.4e, ov_upper: %.4e" %\
                        (global_minimum_viscosity, global_minimum_viscosity, slab_core_viscosity, global_minimum_viscosity, global_minimum_viscosity)
        # 4. The peierls rheology
        if  if_peierls:
            o_dict['Material model'][material_model_subsection]['Include Peierls creep'] = 'true'
            if fix_peierls_V_as is not "":
                o_dict['Material model']['Visco Plastic TwoD']['Activation volume differences for Peierls creep'] =\
                    '%.4e' % rheology[fix_peierls_V_as + '_creep']['V']
        # rewrite the reset viscosity part
        if _type == 's07':
            o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
            "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CV=1e20, Wp=%.4e" % (trailing_length, box_height, box_length, sp_width)
        elif _type == 's07_newton' and reset_trailing_ov_viscosity == 1:
            o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
            "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CV=1e20, Wp=%.4e" % (trailing_length, box_height, box_length, sp_width)
        elif _type == '2d_consistent':
            # selection between two methods
            if trailing_length < 1e-6 and trailing_length_1 > 1e-6:
                reset_length = trailing_length_1
            else:
                reset_length = trailing_length
                # another implementation, set trailing_length = 0.0 to turn off the older implementation
            # different geometry
            if geometry == "box":
                o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
            "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CV=1e20, Wp=%.4e" % (reset_length, box_height, box_length, sp_width)
            elif geometry == "chunk":
                o_dict['Material model'][material_model_subsection]['Reset viscosity function']["Coordinate system"] = "spherical"
                o_dict['Material model'][material_model_subsection]['Reset viscosity function']["Variable names"] = "r, phi, theta"
                o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
                    "Depth=1.45e5, Ro=6.371e6, Width=%.4e, PHIM=%.4e, CV=1e20, Wp=%.4e, PI=3.1415" % (reset_length, box_length/Ro, sp_width)
                o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function expression'] =\
                    "(((r > Ro - Depth) && ((Ro*phi < Width)||(Ro*(PHIM-phi) < Width)) && (Ro*(PI/2.0 - theta) <= Wp))? CV: -1.0)"
        # rewrite the reset density part
        if _type == '2d_consistent':
            if "Reset density function" in o_dict['Material model'][material_model_subsection]:
                # selection between two method
                if trailing_length < 1e-6 and trailing_length_1 > 1e-6:
                    reset_length = trailing_length_1
                else:
                    reset_length = trailing_length
                if geometry == "box":
                    o_dict['Material model'][material_model_subsection]["Reset density function"]['Function constants'] =\
                        "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CD=3300.0, Wp=%.4e" % (reset_length, box_height, box_length, sp_width)
                elif geometry == "chunk":
                    o_dict['Material model'][material_model_subsection]['Reset density function']["Coordinate system"] = "spherical"
                    o_dict['Material model'][material_model_subsection]['Reset density function']["Variable names"] = "r, phi, theta"
                    o_dict['Material model'][material_model_subsection]['Reset density function']['Function constants'] =\
                        "Depth=1.45e5, Ro=6.371e6, Width=%.4e, PHIM=%.4e, CD=3300.0, Wp=%.4e, PI=3.1415" % (reset_length, box_length/Ro, sp_width)
                    o_dict['Material model'][material_model_subsection]['Reset density function']['Function expression'] =\
                        "(((r > Ro - Depth) && ((Ro*phi < Width)||(Ro*(PHIM-phi) < Width)) && (Ro*(PI/2.0 - theta) <= Wp))? CD: -1.0)"

        # rewrite the reaction morb part
        if _type in ['s07', 's07_newton', '2d_consistent']:
            if abs(Dsz - 50e3)/50e3 > 1e-6:
                # for the sake of backwards compatible
                Dsz_str = str(Dsz)
            else:
                Dsz_str = "5e+04"
            if abs(Dsz + Ddl - 100e3)/100e3 > 1e-6:
                # for the sake of backwards compatible
                Dp_str = str(Dsz + Ddl)
            else:
                Dp_str = "1e5"
        if reset_trailing_morb == 1:
            o_dict['Material model'][material_model_subsection]['Reaction mor'] = 'true'
        else:
            o_dict['Material model'][material_model_subsection]['Reaction mor'] = 'false'
        if _type == "s07":
            # fix the variables type-wise
            if reset_trailing_morb == 1:
                o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                    "Width=%.4e, Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5" %  (trailing_length, box_height, box_length, Dsz_str, Dp_str, sp_width)
        elif _type == "s07_newton":
            if reset_trailing_morb == 1:
                o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                    "Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5" %  (box_height, box_length, Dsz_str, Dp_str, sp_width)
        elif _type == '2d_consistent':
            if reset_trailing_morb == 1:
                if geometry == "box":
                    if sp_ridge_x > 0.0 and ov_side_dist > 0.0:
                        o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                            "Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5, Dplate=200e3, Wweak=55e3, pRidge=%.4e, dOvSide=%.4e" % \
                            (box_height, box_length, Dsz_str, Dp_str, sp_width, sp_ridge_x, ov_side_dist)
                elif geometry == "chunk":
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']["Coordinate system"] = "spherical"
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']["Variable names"] = "r, phi, theta"
                    if sp_ridge_x > 0.0 and ov_side_dist > 0.0:
                        o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                            "Ro=6.371e6, PHIM=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5, Dplate=200e3, Wweak=55e3, pRidge=%.4e, dOvSide=%.4e" % \
                            (box_length/Ro, Dsz_str, Dp_str, sp_width, sp_ridge_x, ov_side_dist)
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function expression'] =\
                        "((r > Ro - DpUp) && (Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth) && (Ro*(PI/2.0 - theta) <= Wp)) ? 0:\\\n\
                                        (((r <= Ro - DpUp) && (r > Ro - Dp) && (Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth) && (Ro*(PI/2.0 - theta) <= Wp)) ? 1:\\\n\
                                        ((r > Ro - Dplate) && (((Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth)) || ((Ro*(PHIM-phi) < dOvSide + pWidth) && (Ro*(PHIM-phi) >= dOvSide))) && (Ro*(PI/2.0 - theta) > Wp) && (Ro*(PI/2.0 - theta) <= Wp + Wweak))? 2:-1)"
            elif reset_trailing_morb == 2:
                # use this option to reset the morb composition when sp ridge is 0.0
                if geometry == "box":
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                        "Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5, Dplate=200e3, Wweak=55e3, pRidge=%.4e, dOvSide=%.4e" % \
                        (box_height, box_length, Dsz_str, Dp_str, sp_width, 6e5, 6e5)
                elif geometry == "chunk":
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']["Coordinate system"] = "spherical"
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']["Variable names"] = "r, phi, theta"
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                        "Ro=6.371e6, PHIM=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5, Dplate=200e3, Wweak=55e3, pRidge=%.4e, dOvSide=%.4e" % \
                        (box_length/Ro, Dsz_str, Dp_str, sp_width, 6e5, 6e5)
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function expression'] =\
                        "((r > Ro - DpUp) && (Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth) && (Ro*(PI/2.0 - theta) <= Wp)) ? 0:\\\n\
                                        (((r <= Ro - DpUp) && (r > Ro - Dp) && (Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth) && (Ro*(PI/2.0 - theta) <= Wp)) ? 1:\\\n\
                                        ((r > Ro - Dplate) && (((Ro*phi > pRidge) && (Ro*phi <= pRidge + pWidth)) || ((Ro*(PHIM-phi) < dOvSide + pWidth) && (Ro*(PHIM-phi) >= dOvSide))) && (Ro*(PI/2.0 - theta) > Wp) && (Ro*(PI/2.0 - theta) <= Wp + Wweak))? 2:-1)"
            elif reset_trailing_morb == 3:
                o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                "Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5, Dplate=200e3, Wweak=55e3, pRidge=6.0000e+05, dOvSide=6.0000e+05" %\
                (box_height, box_length, Dsz_str, Dp_str, sp_width)
        else:
            raise ValueError()
        # prescribe mantle temperature
        if _type == '2d_consistent':
            if geometry == "chunk":
                # modify for the chunk geometry
                o_dict["Prescribed temperatures"]["Indicator function"]["Coordinate system"] = "spherical"
                o_dict["Prescribed temperatures"]["Indicator function"]["Variable names"] = "r, phi, theta"
                o_dict["Prescribed temperatures"]["Indicator function"]["Function constants"] = "Depth=1.45e5, Width=9.0000e+05, Wp=1.0000e+06, Ro=6.371e6, PHIM=%.4f, PI=3.1415" % (box_length/Ro)
                o_dict["Prescribed temperatures"]["Indicator function"]["Function expression"] = "(((r>Ro-Depth)&&(Ro*(PI/2.0 - theta)<Wp)&&((Ro*phi<Width)||(Ro*(PHIM-phi)<Width))) ? 1:0)"
            # only do this for the 2d_consistent model
            if prescribe_mantle_sp and prescribe_mantle_ov:
                o_dict["Prescribe internal mantle adiabat temperatures"] = 'true'
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function constants"] = \
                    "Do=%.4e, xm=%.4e, Wp=%.4e, Depth=150e3, pRidge=%.4e, dOvSide=%.4e" % \
                    (box_height, box_length, sp_width, sp_ridge_x, ov_side_dist)
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function expression"] = \
                    "(z > Do - Depth) && ((x <= pRidge) || (x > (xm - dOvSide))) && (y <= Wp)? 1.0: 0.0"
            elif prescribe_mantle_sp and (not prescribe_mantle_ov):
                o_dict["Prescribe internal mantle adiabat temperatures"] = 'true'
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function constants"] = \
                    "Do=%.4e, xm=%.4e, Wp=%.4e, Depth=150e3, pRidge=%.4e" % \
                    (box_height, box_length, sp_width, sp_ridge_x)
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function expression"] = \
                    "(z > Do - Depth) && (x <= pRidge) && (y <= Wp)? 1.0: 0.0"
                pass
            elif (not prescribe_mantle_sp) and prescribe_mantle_ov:
                o_dict["Prescribe internal mantle adiabat temperatures"] = 'true'
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function constants"] = \
                    "Do=%.4e, xm=%.4e, Wp=%.4e, Depth=150e3, dOvSide=%.4e" % \
                    (box_height, box_length, sp_width, ov_side_dist)
                o_dict["Prescribed mantle adiabat temperatures"]["Indicator function"]["Function expression"] = \
                    "(z > Do - Depth) && (x > (xm - dOvSide)) && (y <= Wp)? 1.0: 0.0"
                pass
            else:
                # the default in the file should be false if there is any
                pass
            # modify the prescribe temperature feature if another implementation of plate is applied
            if geometry == "box":
                if trailing_length < 1e-6 and trailing_length_1 > 1e-6:
                    # another implementation, set trailing_length = 0.0 to turn off the older implementation
                    # here we set this option to false and remove the corresponding section
                    Do_str = "2.890e6"
                    if abs(box_height - 2.89e6) / 2.89e6 > 1e-6:
                        Do_str = "%.4e" % box_height
                    prescribe_T_area_width = trailing_length_1 + 300e3
                    o_dict["Prescribe internal mantle adiabat temperatures"] = 'false'
                    if "Prescribed mantle adiabat temperatures" in o_dict:
                        o_dict.pop("Prescribed mantle adiabat temperatures")
                    o_dict["Prescribe internal temperatures"] = 'true'
                    o_dict["Prescribed temperatures"]["Indicator function"]["Function constants"] = \
                        "Depth=1.45e5, Width=%.4e, Do=%s, xm=%.4e, Wp=%.4e" % (prescribe_T_area_width, Do_str, box_length, sp_width)
                    o_dict["Prescribed temperatures"]["Plate model 1"] = {
                        'Area width': "%.4e" % prescribe_T_area_width,
                        'Subducting plate velocity': "%.4e" % (sp_rate / year),
                        'Overiding plate age': "%.4e" % (ov_age * year),
                        "Overiding area width": "%.4e" % prescribe_T_area_width,
                        "Top temperature": "273.0"
                    }
            elif geometry == "chunk":
                if trailing_length < 1e-6 and trailing_length_1 > 1e-6:
                    # another implementation, set trailing_length = 0.0 to turn off the older implementation
                    # here we set this option to false and remove the corresponding section
                    prescribe_T_area_width = trailing_length_1 + 300e3
                    o_dict["Prescribe internal mantle adiabat temperatures"] = 'false'
                    if "Prescribed mantle adiabat temperatures" in o_dict:
                        o_dict.pop("Prescribed mantle adiabat temperatures")
                    o_dict["Prescribe internal temperatures"] = 'true'
                    o_dict["Prescribed temperatures"]["Plate model 1"] = {
                        'Area width': "%.4e" % prescribe_T_area_width,
                        'Subducting plate velocity': "%.4e" % (sp_rate / year),
                        'Overiding plate age': "%.4e" % (ov_age * year),
                        "Overiding area width": "%.4e" % prescribe_T_area_width,
                        "Top temperature": "273.0"
                    }
        # reset composition viscosity
        # usage is to prescribe a strong core with the plate
        if reset_composition_viscosity:
            o_dict['Material model'][material_model_subsection]['Reset composition viscosity'] = 'true'
            function_dict = o_dict['Material model'][material_model_subsection]['Reset composition viscosity function']
            if _type == "s07":
                function_dict["Function constants"] = "Width=%.4e, CV=1e23" % reset_composition_viscosity_width
            else:
                raise NotImplementedError()
            o_dict['Material model'][material_model_subsection]['Reset composition viscosity function'] = function_dict



        o_dict['Material model'][material_model_subsection] = {**o_dict['Material model'][material_model_subsection], **outputs}  # prepare entries

        # solver options
        if stokes_solver_type == "block AMG":
            # here the default value is AMG, if not present, don't change anything
            try:
                if "Stokes solver type" in o_dict['Solver parameters']["Stokes solver parameters"]:
                    o_dict['Solver parameters']['Stokes solver parameters']["Stokes solver type"] = "block AMG"
                if "Skip expensive stokes solver" in o_dict['Solver parameters']["Stokes solver parameters"]:
                    o_dict['Solver parameters']["Stokes solver parameters"].pop("Skip expensive stokes solver")
            except KeyError:
                pass
        elif stokes_solver_type == "block GMG":
            o_dict['Solver parameters']['Stokes solver parameters']["Stokes solver type"] = "block GMG"
        elif stokes_solver_type == "block GMG with iterated defect correction Stokes":
            o_dict['Nonlinear solver scheme'] = "single Advection, iterated defect correction Stokes"
            o_dict['Solver parameters']['Stokes solver parameters']["Stokes solver type"] = "block GMG"
        else:
            raise ValueError("stokes_solver_type must be in [block AMG, block GMG].")

        self.idict = o_dict
        pass


    def configure_wb(self, _type, if_wb, geometry, sp_width, sp_length, trailing_length, Dsz, Ddl, slab_length,\
    dip_angle, sp_age_trench, ov_age, setup_method, sp_rate, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, make_2d_consistent_plate):
        '''
        Configure wb file
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # header options
        if geometry == 'chunk':
            self.wb_dict["coordinate system"] = \
            {
                "model": "spherical",
                "depth method": "begin at end segment"
            }
        # geometry options
        if setup_method == 'manual':
            if _type in ["s07", "s07T"]:
                wb_configure_plate_schellart07(self.wb_dict, sp_width, sp_length, trailing_length, Dsz, Ddl, slab_length)
            elif _type in ["s07_newton"]:
                wb_configure_plate_schellart07_Tdependent(self.wb_dict, sp_width, sp_length, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age)
            else:
                raise ValueError("Wrong value for \"type\"")
        elif setup_method == '2d_consistent':
            # first select between different geometries
            if geometry == 'chunk':
                if make_2d_consistent_plate == 2:
                    self.wb_dict = wb_configure_plate_2d_consistent_2_chunk(self.wb_dict, sp_width, sp_rate, Dsz,\
                        slab_length, dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate,  if_ov_trans,\
                        ov_trans_age, ov_trans_length, sp_ridge_x, ov_side_dist, box_length)
                else:
                    raise ValueError("chunk geometry has implemented this option of make_2d_consistent_plate %d" % make_2d_consistent_plate)
            elif geometry == 'box':
                if make_2d_consistent_plate == 1:
                    self.wb_dict = wb_configure_plate_2d_consistent_1(self.wb_dict, sp_width, sp_rate, Dsz,\
                        slab_length, dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate,  if_ov_trans,\
                        ov_trans_age, ov_trans_length, sp_ridge_x, ov_side_dist, box_length)
                elif make_2d_consistent_plate == 2:
                    self.wb_dict = wb_configure_plate_2d_consistent_2(self.wb_dict, sp_width, sp_rate, Dsz,\
                        slab_length, dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate,  if_ov_trans,\
                        ov_trans_age, ov_trans_length, sp_ridge_x, ov_side_dist, box_length)
                elif make_2d_consistent_plate == 3:
                    self.wb_dict = wb_configure_plate_2d_consistent_3(self.wb_dict, sp_width, sp_rate, Dsz,\
                        slab_length, dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate,  if_ov_trans,\
                        ov_trans_age, ov_trans_length, sp_ridge_x, ov_side_dist, box_length)
                else:
                    self.wb_dict = wb_configure_plate_2d_consistent(self.wb_dict, sp_width, sp_rate, Dsz, Ddl,\
                        slab_length, dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate,  if_ov_trans,\
                        ov_trans_age, ov_trans_length, sp_ridge_x, ov_side_dist, box_length)
            else:
                raise ValueError("Geometry must by \"chunk\" or \"box\", get %s" % geometry)
            pass
    

def wb_configure_plate_schellart07(wb_dict, sp_width, sp_length, trailing_width, Dsz, Ddl, slab_length):
    '''
    World builder configuration of plates in Schellart etal 2007
    '''
    o_dict = wb_dict.copy()
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[trailing_width, -sp_width], [trailing_width, sp_width], [sp_length, sp_width] ,[sp_length, -sp_width]]
    if abs(Dsz - 50e3) / 50e3 > 1e-6:
        sp_dict["composition models"][0]["max depth"] = Dsz
        sp_dict["composition models"][1]["min depth"] = Dsz
    if abs(Dsz + Ddl - 100e3) / 100e3 > 1e-6:
        sp_dict["composition models"][1]["max depth"] = Dsz + Ddl
    o_dict['features'][i0] = sp_dict
    # slab
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    sdict = o_dict['features'][i0]
    sdict["coordinates"] = [[sp_length, -sp_width], [sp_length, sp_width]]
    if abs(Dsz - 50e3) / 50e3 > 1e-6:
        for i in range(len(sdict["segments"])-1):
            # the last one is a phantom for temperature tapering
            segment = sdict["segments"][i]
            segment["composition models"][0]["max distance slab top"] = Dsz
            segment["composition models"][1]["min distance slab top"] = Dsz
            sdict["segments"][i] = segment
    if abs(Dsz + Ddl - 100e3) / 100e3 > 1e-6:
        for i in range(len(sdict["segments"])-1):
            segment = sdict["segments"][i]
            segment["composition models"][1]["max distance slab top"] = Dsz + Ddl
            sdict["segments"][i] = segment
    if abs(slab_length - 167e3)/167e3 > 1e-6:
        segment = sdict["segments"][1]
        segment['length'] = slab_length
        sdict["segments"][1] = segment
    o_dict['features'][i0] = sdict


def wb_configure_plate_schellart07_Tdependent(wb_dict, sp_width, sp_length, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age):
    '''
    World builder configuration of plates in Schellart etal 2007
    '''
    o_dict = wb_dict.copy()
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if abs(ov_age - 40e6) / 40e6 > 1e-6:
        ov_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ov_dict
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[0.0, -sp_width], [0.0, sp_width], [sp_length, sp_width] ,[sp_length, -sp_width]]
    if abs(Dsz - 50e3) / 50e3 > 1e-6:
        sp_dict["composition models"][0]["max depth"] = Dsz
        sp_dict["composition models"][1]["min depth"] = Dsz
    if abs(Dsz + Ddl - 100e3) / 100e3 > 1e-6:
        sp_dict["composition models"][1]["max depth"] = Dsz + Ddl
    o_dict['features'][i0] = sp_dict
    # slab
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    sdict = o_dict['features'][i0]
    sdict["coordinates"] = [[sp_length, -sp_width], [sp_length, sp_width]]
    if abs(Dsz - 50e3) / 50e3 > 1e-6:  # moho
        for i in range(len(sdict["segments"])-1):
            # the last one is a phantom for temperature tapering
            segment = sdict["segments"][i]
            segment["composition models"][0]["max distance slab top"] = Dsz
            segment["composition models"][1]["min distance slab top"] = Dsz
            sdict["segments"][i] = segment
    if abs(Dsz + Ddl - 100e3) / 100e3 > 1e-6: # lithosphere - athenosphere boundary
        for i in range(len(sdict["segments"])-1):
            segment = sdict["segments"][i]
            segment["composition models"][1]["max distance slab top"] = Dsz + Ddl
            sdict["segments"][i] = segment
    if abs(slab_length - 167e3)/167e3 > 1e-6:
        segment = sdict["segments"][1]
        segment['length'] = slab_length
        sdict["segments"][1] = segment
    if abs(dip_angle - 15.5)/15.5 > 1e-6:
        # initial dipping angle
        segment = sdict["segments"][0]
        segment['angle'] = [0, dip_angle]
        segment['length'] = 20e3 * dip_angle/15.5 # maintain the same curvature
        segment = sdict["segments"][1]
        segment['angle'] = [dip_angle, dip_angle]
        segment = sdict["segments"][2]
        segment['angle'] = [dip_angle, dip_angle]
    o_dict['features'][i0] = sdict


def wb_configure_plate_2d_consistent(wb_dict, sp_width, sp_rate, Dsz, Ddl, slab_length,\
    dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, **kwargs):
    '''
    World builder configuration of plates in Schellart etal 2007
    '''
    Xmax = 40000e3  # some big values that are not ever reached
    Ymax = 20000e3
    pe_width = 55e3 # width of the plate edges
    Ro = kwargs.get("Ro", 6371e3)
    sp_length = sp_rate * sp_age_trench
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [Xmax, 0.0]]
    # overiding plate 1
    if if_ov_trans and ov_age > (1e6 + ov_trans_age):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_feature, ov =\
            wb_configure_transit_ov_plates(wb_dict['features'][i0], sp_length + sp_ridge_x,\
                ov_age, ov_trans_age, ov_trans_length, wb_new_ridge,\
                Ro=Ro, geometry='box')
        o_dict['features'][i0] = ov_trans_feature
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except ParsePrm.WBFeatureNotFoundError:
            pass
        else:
            o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = sp_length + sp_ridge_x
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if ov_side_dist > 0.0:
        ov_end = box_length - ov_side_dist
    else:
        ov_end = Xmax
    ov_dict["coordinates"] = [[ov, -sp_width], [ov, sp_width], [ov_end, sp_width] ,[ov_end, -sp_width]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_ridge_x + sp_length, sp_width], [sp_ridge_x + sp_length, sp_width + pe_width], [ov_end, sp_width + pe_width] ,[ov_end, sp_width]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[sp_ridge_x, -sp_width], [sp_ridge_x, sp_width], [sp_ridge_x + sp_length, sp_width] ,[sp_ridge_x + sp_length, -sp_width]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz + Ddl  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[sp_ridge_x, sp_width], [sp_ridge_x, sp_width + pe_width], [sp_ridge_x + sp_length, sp_width + pe_width] ,[sp_ridge_x + sp_length, sp_width]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    spe_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    # side plate
    if assign_side_plate == 1:
        sdp_dict = deepcopy(ov_dict)
        sdp_dict["name"] = "side plate"
        sdp_dict["coordinates"] = [[0.0, sp_width+pe_width], [0.0, Ymax], [Xmax, Ymax] ,[Xmax, sp_width+pe_width]]
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, "side plate")
        except ParsePrm.WBFeatureNotFoundError:
            i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
            o_dict['features'].insert(i0 + 1, sdp_dict)
        else:
            o_dict['features'][i0] = sdp_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length + sp_ridge_x, -sp_width], [sp_length + sp_ridge_x, sp_width]]
    sdict["dip point"] = [Xmax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz + Ddl
        sdict["segments"][i] = segment
    if wb_new_ridge == 1:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]]
    else:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]
    sdict["temperature models"][0]["plate velocity"] = sp_rate
    return o_dict


def wb_configure_plate_2d_consistent_1(wb_dict, sp_width, sp_rate, Dsz, slab_length,\
    dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, **kwargs):
    '''
    setup up consistently with the 2d models:
    a. scale the harzburgite layer with the crustal layer.
    '''
    Xmax = 40000e3  # some big values that are not ever reached
    Ymax = 20000e3
    D2C_ratio = 35.2e3 / 7.5e3 # ratio of depleted / crust layer
    pe_width = 55e3 # width of the plate edges
    Ro = kwargs.get("Ro", 6371e3)
    sp_length = sp_rate * sp_age_trench
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [Xmax, 0.0]]
    # overiding plate 1
    if if_ov_trans and ov_age > (1e6 + ov_trans_age):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_feature, ov =\
            wb_configure_transit_ov_plates(wb_dict['features'][i0], sp_length + sp_ridge_x,\
                ov_age, ov_trans_age, ov_trans_length, wb_new_ridge,\
                Ro=Ro, geometry='box')
        o_dict['features'][i0] = ov_trans_feature
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except ParsePrm.WBFeatureNotFoundError:
            pass
        else:
            o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = sp_length + sp_ridge_x
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if ov_side_dist > 0.0:
        ov_end = box_length - ov_side_dist
    else:
        ov_end = Xmax
    ov_dict["coordinates"] = [[ov, -sp_width], [ov, sp_width], [ov_end, sp_width] ,[ov_end, -sp_width]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_ridge_x + sp_length, sp_width], [sp_ridge_x + sp_length, sp_width + pe_width], [ov_end, sp_width + pe_width] ,[ov_end, sp_width]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[sp_ridge_x, -sp_width], [sp_ridge_x, sp_width], [sp_ridge_x + sp_length, sp_width] ,[sp_ridge_x + sp_length, -sp_width]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[sp_ridge_x, sp_width], [sp_ridge_x, sp_width + pe_width], [sp_ridge_x + sp_length, sp_width + pe_width] ,[sp_ridge_x + sp_length, sp_width]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    spe_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    # side plate
    if assign_side_plate == 1:
        sdp_dict = deepcopy(ov_dict)
        sdp_dict["name"] = "side plate"
        sdp_dict["coordinates"] = [[0.0, sp_width+pe_width], [0.0, Ymax], [Xmax, Ymax] ,[Xmax, sp_width+pe_width]]
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, "side plate")
        except ParsePrm.WBFeatureNotFoundError:
            i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
            o_dict['features'].insert(i0 + 1, sdp_dict)
        else:
            o_dict['features'][i0] = sdp_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length + sp_ridge_x, -sp_width], [sp_length + sp_ridge_x, sp_width]]
    sdict["dip point"] = [Xmax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz * D2C_ratio
        sdict["segments"][i] = segment
    if wb_new_ridge == 1:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]]
    else:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]
    sdict["temperature models"][0]["plate velocity"] = sp_rate
    return o_dict

def wb_configure_plate_2d_consistent_2(wb_dict, sp_width, sp_rate, Dsz, slab_length,\
    dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, **kwargs):
    '''
    setup up consistently with the 2d models:
    a. scale the harzburgite layer with the crustal layer.
    version 3: add the composition of the overriding plate
    '''
    Xmax = 40000e3  # some big values that are not ever reached
    Ymax = 20000e3
    D2C_ratio = 35.2e3 / 7.5e3 # ratio of depleted / crust layer
    pe_width = 55e3 # width of the plate edges
    Ro = kwargs.get("Ro", 6371e3)
    sp_length = sp_rate * sp_age_trench
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [Xmax, 0.0]]
    # overiding plate 1
    if if_ov_trans and ov_age > (1e6 + ov_trans_age):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_dict, ov =\
            wb_configure_transit_ov_plates_1(wb_dict['features'][i0], sp_length + sp_ridge_x,\
                ov_age, ov_trans_age, ov_trans_length, wb_new_ridge, sp_width,\
                Ro=Ro, geometry='box')
        ov_trans_dict["composition models"][0]["max depth"] = Dsz  # moho
        ov_trans_dict["composition models"][1]["min depth"] = Dsz
        ov_trans_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
        o_dict['features'][i0] = ov_trans_dict
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except ParsePrm.WBFeatureNotFoundError:
            pass
        else:
            o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = sp_length + sp_ridge_x
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if ov_side_dist > 0.0:
        ov_end = box_length - ov_side_dist
    else:
        ov_end = Xmax
    ov_dict["coordinates"] = [[ov, -sp_width], [ov, sp_width], [ov_end, sp_width] ,[ov_end, -sp_width]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    ov_dict["composition models"][0]["max depth"] = Dsz  # moho
    ov_dict["composition models"][1]["min depth"] = Dsz
    ov_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_ridge_x + sp_length, sp_width], [sp_ridge_x + sp_length, sp_width + pe_width], [ov_end, sp_width + pe_width] ,[ov_end, sp_width]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[sp_ridge_x, -sp_width], [sp_ridge_x, sp_width], [sp_ridge_x + sp_length, sp_width] ,[sp_ridge_x + sp_length, -sp_width]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[sp_ridge_x, sp_width], [sp_ridge_x, sp_width + pe_width], [sp_ridge_x + sp_length, sp_width + pe_width] ,[sp_ridge_x + sp_length, sp_width]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    spe_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    # side plate
    if assign_side_plate == 1:
        sdp_dict = deepcopy(ov_dict)
        sdp_dict["name"] = "side plate"
        sdp_dict["coordinates"] = [[0.0, sp_width+pe_width], [0.0, Ymax], [Xmax, Ymax] ,[Xmax, sp_width+pe_width]]
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, "side plate")
        except ParsePrm.WBFeatureNotFoundError:
            i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
            o_dict['features'].insert(i0 + 1, sdp_dict)
        else:
            o_dict['features'][i0] = sdp_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length + sp_ridge_x, -sp_width], [sp_length + sp_ridge_x, sp_width]]
    sdict["dip point"] = [Xmax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz * D2C_ratio
        sdict["segments"][i] = segment
    if wb_new_ridge == 1:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]]
    else:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]
    sdict["temperature models"][0]["plate velocity"] = sp_rate
    return o_dict


def wb_configure_plate_2d_consistent_2_chunk(wb_dict, sp_width, sp_rate, Dsz, slab_length,\
    dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, **kwargs):
    '''
    setup up consistently with the 2d models:
    a. scale the harzburgite layer with the crustal layer.
    version 3: add the composition of the overriding plate
    specific for the chunk geometry
    '''
    Xmax = 40000e3  # some big values that are not ever reached
    PHImax = 360.0
    Ymax = 20000e3
    THETAmax = 180.0
    D2C_ratio = 35.2e3 / 7.5e3 # ratio of depleted / crust layer
    pe_width = 55e3 # width of the plate edges
    Ro = kwargs.get("Ro", 6371e3)
    sp_width_deg = sp_width / Ro * 180.0 / np.pi
    print("sp_width: ", sp_width) # debug
    print("sp_width_deg: ", sp_width_deg)
    sp_ridge_x_deg = sp_ridge_x / Ro * 180.0 / np.pi
    sp_length = sp_rate * sp_age_trench
    sp_length_deg = sp_length / Ro * 180.0 / np.pi
    pe_width_deg = pe_width / Ro * 180.0 / np.pi
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [PHImax, 0.0]]
    # overiding plate 1
    if if_ov_trans and ov_age > (1e6 + ov_trans_age):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_dict, ov =\
            wb_configure_transit_ov_plates_1_chunk(wb_dict['features'][i0], sp_length + sp_ridge_x,\
                ov_age, ov_trans_age, ov_trans_length, wb_new_ridge, sp_width,\
                Ro=Ro, geometry='box')
        ov_trans_dict["composition models"][0]["max depth"] = Dsz  # moho
        ov_trans_dict["composition models"][1]["min depth"] = Dsz
        ov_trans_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
        o_dict['features'][i0] = ov_trans_dict
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except ParsePrm.WBFeatureNotFoundError:
            pass
        else:
            o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = sp_length + sp_ridge_x
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if ov_side_dist > 0.0:
        ov_end = box_length - ov_side_dist
        ov_end_deg = ov_end / Ro * 180.0 / np.pi
    else:
        ov_end = Xmax
        ov_end_deg = 360.0
    ov_deg = ov / Ro * 180.0 / np.pi
    ov_dict["coordinates"] = [[ov_deg, -1.0], [ov_deg, sp_width_deg], [ov_end_deg, sp_width_deg] ,[ov_end_deg, -1.0]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    ov_dict["composition models"][0]["max depth"] = Dsz  # moho
    ov_dict["composition models"][1]["min depth"] = Dsz
    ov_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_ridge_x_deg + sp_length_deg, sp_width_deg], [sp_ridge_x_deg + sp_length_deg, sp_width_deg + pe_width_deg], [ov_end_deg, sp_width_deg + pe_width_deg] ,[ov_end_deg, sp_width_deg]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[sp_ridge_x_deg, -1.0], [sp_ridge_x_deg, sp_width_deg], [sp_ridge_x_deg + sp_length_deg, sp_width_deg] ,[sp_ridge_x_deg + sp_length_deg, -1.0]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x_deg, -1.0], [sp_ridge_x_deg, 180.0]]]
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[sp_ridge_x_deg, sp_width_deg], [sp_ridge_x_deg, sp_width_deg + pe_width_deg], [sp_ridge_x_deg + sp_length_deg, sp_width_deg + pe_width_deg] ,[sp_ridge_x_deg + sp_length_deg, sp_width_deg]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    spe_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x_deg, -1.0], [sp_ridge_x_deg, 180.0]]]
    # side plate
    if assign_side_plate == 1:
        sdp_dict = deepcopy(ov_dict)
        sdp_dict["name"] = "side plate"
        sdp_dict["coordinates"] = [[0.0, sp_width_deg+pe_width_deg], [0.0, THETAmax], [PHImax, THETAmax] ,[PHImax, sp_width_deg+pe_width_deg]]
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, "side plate")
        except ParsePrm.WBFeatureNotFoundError:
            i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
            o_dict['features'].insert(i0 + 1, sdp_dict)
        else:
            o_dict['features'][i0] = sdp_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length_deg + sp_ridge_x_deg, -1.0], [sp_length_deg + sp_ridge_x_deg, sp_width_deg]]
    sdict["dip point"] = [PHImax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz * D2C_ratio
        sdict["segments"][i] = segment
    if wb_new_ridge == 1:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[[sp_ridge_x_deg,-1.0], [sp_ridge_x_deg, 180.0]]]
    else:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[sp_ridge_x_deg,-1.0], [sp_ridge_x_deg, 180.0]]
    sdict["temperature models"][0]["plate velocity"] = sp_rate
    # mantle to subtract
    i0 = ParsePrm.FindWBFeatures(o_dict, "mantle to substract")
    mdict = o_dict['features'][i0] 
    mdict["coordinates"] = [[0.0, -1.0], [0.0, THETAmax], [PHImax,THETAmax], [PHImax, -1.0]]
    
    return o_dict

# todo_ov
def wb_configure_plate_2d_consistent_3(wb_dict, sp_width, sp_rate, Dsz, slab_length,\
    dip_angle, sp_age_trench, ov_age, wb_new_ridge, assign_side_plate, if_ov_trans, ov_trans_age,\
    ov_trans_length, sp_ridge_x, ov_side_dist, box_length, **kwargs):
    '''
    setup up consistently with the 2d models:
    a. scale the harzburgite layer with the crustal layer.
    version 3: remove the composition of the overriding plate
    '''
    Xmax = 40000e3  # some big values that are not ever reached
    Ymax = 20000e3
    D2C_ratio = 35.2e3 / 7.5e3 # ratio of depleted / crust layer
    pe_width = 55e3 # width of the plate edges
    Ro = kwargs.get("Ro", 6371e3)
    sp_length = sp_rate * sp_age_trench
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [Xmax, 0.0]]
    # overiding plate 1
    if if_ov_trans and ov_age > (1e6 + ov_trans_age):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_dict, ov =\
            wb_configure_transit_ov_plates_1(wb_dict['features'][i0], sp_length + sp_ridge_x,\
                ov_age, ov_trans_age, ov_trans_length, wb_new_ridge, sp_width,\
                Ro=Ro, geometry='box')
        if "composition models" in ov_trans_dict:
            ov_trans_dict.pop("composition models")
        o_dict['features'][i0] = ov_trans_dict
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except ParsePrm.WBFeatureNotFoundError:
            pass
        else:
            o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = sp_length + sp_ridge_x
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    if ov_side_dist > 0.0:
        ov_end = box_length - ov_side_dist
    else:
        ov_end = Xmax
    ov_dict["coordinates"] = [[ov, -sp_width], [ov, sp_width], [ov_end, sp_width] ,[ov_end, -sp_width]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    if "composition models" in ov_dict:
        ov_dict.pop("composition models")
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_ridge_x + sp_length, sp_width], [sp_ridge_x + sp_length, sp_width + pe_width], [ov_end, sp_width + pe_width] ,[ov_end, sp_width]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[sp_ridge_x, -sp_width], [sp_ridge_x, sp_width], [sp_ridge_x + sp_length, sp_width] ,[sp_ridge_x + sp_length, -sp_width]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz * D2C_ratio  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[sp_ridge_x, sp_width], [sp_ridge_x, sp_width + pe_width], [sp_ridge_x + sp_length, sp_width + pe_width] ,[sp_ridge_x + sp_length, sp_width]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    spe_dict["temperature models"][0]["ridge coordinates"] = [[[sp_ridge_x, -10000000.0], [sp_ridge_x, 10000000.0]]]
    # side plate
    if assign_side_plate == 1:
        sdp_dict = deepcopy(ov_dict)
        if "composition models" in sdp_dict:
            # remove composition
            sdp_dict.pop("composition models")
        sdp_dict["name"] = "side plate"
        sdp_dict["coordinates"] = [[0.0, sp_width+pe_width], [0.0, Ymax], [Xmax, Ymax] ,[Xmax, sp_width+pe_width]]
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, "side plate")
        except ParsePrm.WBFeatureNotFoundError:
            i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
            o_dict['features'].insert(i0 + 1, sdp_dict)
        else:
            o_dict['features'][i0] = sdp_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length + sp_ridge_x, -sp_width], [sp_length + sp_ridge_x, sp_width]]
    sdict["dip point"] = [Xmax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz * D2C_ratio
        sdict["segments"][i] = segment
    if wb_new_ridge == 1:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]]
    else:
        sdict["temperature models"][0]["ridge coordinates"] = \
            [[sp_ridge_x,-10000000.0], [sp_ridge_x, 10000000.0]]
    sdict["temperature models"][0]["plate velocity"] = sp_rate
    return o_dict


def s07T_assign_mantle_rheology(o_dict, rheology):
    '''
    Assign mantle rheology in the s07T model
    Note this model has 2 compositions:
        sp_upper and sp_lower
    ''' 
    diff_crust_A = 5e-21
    diff_crust_m = 0.0
    diff_crust_E = 0.0
    diff_crust_V = 0.0
    disl_crust_A = 5e-32
    disl_crust_n = 1.0
    disl_crust_E = 0.0
    disl_crust_V = 0.0
    diffusion_creep = rheology['diffusion_creep']
    dislocation_creep = rheology['dislocation_creep']
    diffusion_creep_lm = rheology['diffusion_lm']
    diff_A = diffusion_creep['A']
    diff_m = diffusion_creep['m']
    diff_n = diffusion_creep['n']
    diff_E = diffusion_creep['E']
    diff_V = diffusion_creep['V']
    diff_d = diffusion_creep['d']
    disl_A = dislocation_creep['A']
    disl_m = dislocation_creep['m']
    disl_n = dislocation_creep['n']
    disl_E = dislocation_creep['E']
    disl_V = dislocation_creep['V']
    disl_d = dislocation_creep['d']
    diff_A_lm = diffusion_creep_lm['A']
    diff_m_lm = diffusion_creep_lm['m']
    diff_n_lm = diffusion_creep_lm['n']
    diff_E_lm = diffusion_creep_lm['E']
    diff_V_lm = diffusion_creep_lm['V']
    diff_d_lm = diffusion_creep_lm['d']
    disl_A_lm = 5.0000e-32
    disl_n_lm = 1.0
    disl_E_lm = 0.0
    disl_V_lm = 0.0
    o_dict['Material model']['Visco Plastic TwoD']['Prefactors for diffusion creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (diff_A, diff_A_lm, diff_crust_A, diff_A, diff_A_lm, diff_A, diff_A_lm, diff_crust_A, diff_A, diff_A_lm)

    o_dict['Material model']['Visco Plastic TwoD']['Grain size exponents for diffusion creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (diff_m, diff_m_lm, diff_crust_m, diff_m, diff_m_lm, diff_m, diff_m_lm, diff_crust_m, diff_m, diff_m_lm)
    
    o_dict['Material model']['Visco Plastic TwoD']['Activation energies for diffusion creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (diff_E, diff_E_lm, diff_crust_E, diff_E, diff_E_lm, diff_E, diff_E_lm, diff_crust_E, diff_E, diff_E_lm)

    o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for diffusion creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (diff_V, diff_V_lm, diff_crust_V, diff_V, diff_V_lm, diff_V, diff_V_lm, diff_crust_V, diff_V, diff_V_lm)
    
    o_dict['Material model']['Visco Plastic TwoD']['Prefactors for dislocation creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (disl_A, disl_A_lm, disl_crust_A, disl_A, disl_A_lm, disl_A, disl_A_lm, disl_crust_A, disl_A, disl_A_lm)

    o_dict['Material model']['Visco Plastic TwoD']['Stress exponents for dislocation creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (disl_n, disl_n_lm, disl_crust_n, disl_n, disl_n_lm, disl_n, disl_n_lm, disl_crust_n, disl_n, disl_n_lm)
    
    o_dict['Material model']['Visco Plastic TwoD']['Activation energies for dislocation creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (disl_E, disl_E_lm, disl_crust_E, disl_E, disl_E_lm, disl_E, disl_E_lm, disl_crust_E, disl_E, disl_E_lm)

    o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for dislocation creep'] = \
        "background: %.4e|%.4e,\
sp_upper: %.4e|%.4e|%.4e,\
sp_lower: %.4e|%.4e,\
plate_edge: %.4e|%.4e|%.4e"\
        % (disl_V, disl_V_lm, disl_crust_V, disl_V, disl_V_lm, disl_V, disl_V_lm, disl_crust_V, disl_V, disl_V_lm)


def consistent_2d_assign_mantle_rheology(o_dict, rheology, **kwargs):
    '''
    Assign mantle rheology in the s07T model
    Note this model has 2 compositions:
        sp_upper and sp_lower
    Inputs:
        kwargs:
            include_ov_upper_plate
    ''' 
    include_ov_upper_plate = kwargs.get("include_ov_upper_plate", 0)
    diff_crust_A = 5e-21
    diff_crust_m = 0.0
    diff_crust_E = 0.0
    diff_crust_V = 0.0
    disl_crust_A = 5e-32
    disl_crust_n = 1.0
    disl_crust_E = 0.0
    disl_crust_V = 0.0
    diffusion_creep = rheology['diffusion_creep']
    dislocation_creep = rheology['dislocation_creep']
    diffusion_creep_lm = rheology['diffusion_lm']
    diff_A = diffusion_creep['A']
    diff_m = diffusion_creep['m']
    diff_n = diffusion_creep['n']
    diff_E = diffusion_creep['E']
    diff_V = diffusion_creep['V']
    diff_d = diffusion_creep['d']
    disl_A = dislocation_creep['A']
    disl_m = dislocation_creep['m']
    disl_n = dislocation_creep['n']
    disl_E = dislocation_creep['E']
    disl_V = dislocation_creep['V']
    disl_d = dislocation_creep['d']
    diff_A_lm = diffusion_creep_lm['A']
    diff_m_lm = diffusion_creep_lm['m']
    diff_n_lm = diffusion_creep_lm['n']
    diff_E_lm = diffusion_creep_lm['E']
    diff_V_lm = diffusion_creep_lm['V']
    diff_d_lm = diffusion_creep_lm['d']
    disl_A_lm = 5.0000e-32
    disl_n_lm = 1.0
    disl_E_lm = 0.0
    disl_V_lm = 0.0
    if include_ov_upper_plate:
        o_dict['Material model']['Visco Plastic TwoD']['Prefactors for diffusion creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (diff_A,diff_A,diff_A,diff_A,diff_A_lm,diff_A_lm,diff_A_lm,diff_A_lm,\
               diff_crust_A, diff_A, diff_A_lm, diff_A_lm,\
                diff_A,diff_A,diff_A,diff_A,diff_A_lm,diff_A_lm,diff_A_lm,diff_A_lm,\
                    diff_crust_A, diff_A, diff_A_lm,\
                    diff_A, diff_A, diff_A_lm, diff_A_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Grain size exponents for diffusion creep'] = \
            "background:  %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (diff_m,diff_m,diff_m,diff_m,diff_m_lm,diff_m_lm,diff_m_lm,diff_m_lm,\
               diff_crust_m, diff_m, diff_m_lm,diff_m_lm,\
                diff_m,diff_m,diff_m,diff_m,diff_m_lm,diff_m_lm,diff_m_lm,diff_m_lm,\
                    diff_crust_m, diff_m, diff_m_lm,\
                    diff_m, diff_m, diff_m_lm,diff_m_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Activation energies for diffusion creep'] = \
            "background:  %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (diff_E,diff_E,diff_E,diff_E,diff_E_lm,diff_E_lm,diff_E_lm,diff_E_lm,\
                diff_crust_E, diff_E, diff_E_lm,diff_E_lm,\
                    diff_E,diff_E,diff_E,diff_E,diff_E_lm,diff_E_lm,diff_E_lm,diff_E_lm,\
                        diff_crust_E, diff_E, diff_E_lm,\
                        diff_E, diff_E, diff_E_lm,diff_E_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for diffusion creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (diff_V,diff_V,diff_V,diff_V,diff_V_lm,diff_V_lm,diff_V_lm,diff_V_lm,\
                diff_crust_V, diff_V, diff_V_lm,diff_V_lm,\
                      diff_V,diff_V,diff_V,diff_V,diff_V_lm,diff_V_lm,diff_V_lm,diff_V_lm,\
                          diff_crust_V, diff_V, diff_V_lm,\
                          diff_V, diff_V, diff_V_lm,diff_V_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Prefactors for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            %(disl_A,disl_A,disl_A,disl_A,disl_A_lm,disl_A_lm,disl_A_lm,disl_A_lm,\
               disl_crust_A, disl_A, disl_A_lm, disl_A_lm,\
                disl_A,disl_A,disl_A,disl_A,disl_A_lm,disl_A_lm,disl_A_lm,disl_A_lm,\
                    disl_crust_A, disl_A, disl_A_lm,\
                    disl_A, disl_A, disl_A_lm, disl_A_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Stress exponents for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (disl_n,disl_n,disl_n,disl_n,disl_n_lm,disl_n_lm,disl_n_lm,disl_n_lm,\
                disl_crust_n, disl_n, disl_n_lm, disl_n_lm,\
                      disl_n,disl_n,disl_n,disl_n,disl_n_lm,disl_n_lm,disl_n_lm,disl_n_lm,\
                          disl_crust_n, disl_n, disl_n_lm,\
                          disl_n, disl_n, disl_n_lm, disl_n_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Activation energies for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (disl_E,disl_E,disl_E,disl_E,disl_E_lm,disl_E_lm,disl_E_lm,disl_E_lm,\
                disl_crust_E, disl_E, disl_E_lm, disl_E_lm,\
                    disl_E,disl_E,disl_E,disl_E,disl_E_lm,disl_E_lm,disl_E_lm,disl_E_lm,\
                        disl_crust_E, disl_E, disl_E_lm,\
                        disl_E, disl_E, disl_E_lm, disl_E_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e,\
    ov_upper: %.4e|%.4e|%.4e|%.4e"\
            % (disl_V,disl_V,disl_V,disl_V,disl_V_lm,disl_V_lm,disl_V_lm,disl_V_lm,\
                disl_crust_V, disl_V, disl_V_lm, disl_V_lm,\
                    disl_V,disl_V,disl_V,disl_V,disl_V_lm,disl_V_lm,disl_V_lm,disl_V_lm,\
                        disl_crust_V, disl_V, disl_V_lm,\
                        disl_V, disl_V, disl_V_lm, disl_V_lm)
    else:
        o_dict['Material model']['Visco Plastic TwoD']['Prefactors for diffusion creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (diff_A,diff_A,diff_A,diff_A,diff_A_lm,diff_A_lm,diff_A_lm,diff_A_lm,\
               diff_crust_A, diff_A, diff_A_lm, diff_A_lm,\
                diff_A,diff_A,diff_A,diff_A,diff_A_lm,diff_A_lm,diff_A_lm,diff_A_lm,\
                    diff_crust_A, diff_A, diff_A_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Grain size exponents for diffusion creep'] = \
            "background:  %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (diff_m,diff_m,diff_m,diff_m,diff_m_lm,diff_m_lm,diff_m_lm,diff_m_lm,\
               diff_crust_m, diff_m, diff_m_lm,diff_m_lm,\
                diff_m,diff_m,diff_m,diff_m,diff_m_lm,diff_m_lm,diff_m_lm,diff_m_lm,\
                    diff_crust_m, diff_m, diff_m_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Activation energies for diffusion creep'] = \
            "background:  %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (diff_E,diff_E,diff_E,diff_E,diff_E_lm,diff_E_lm,diff_E_lm,diff_E_lm,\
                diff_crust_E, diff_E, diff_E_lm,diff_E_lm,\
                    diff_E,diff_E,diff_E,diff_E,diff_E_lm,diff_E_lm,diff_E_lm,diff_E_lm,\
                        diff_crust_E, diff_E, diff_E_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for diffusion creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (diff_V,diff_V,diff_V,diff_V,diff_V_lm,diff_V_lm,diff_V_lm,diff_V_lm,\
                diff_crust_V, diff_V, diff_V_lm,diff_V_lm,\
                      diff_V,diff_V,diff_V,diff_V,diff_V_lm,diff_V_lm,diff_V_lm,diff_V_lm,\
                          diff_crust_V, diff_V, diff_V_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Prefactors for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            %(disl_A,disl_A,disl_A,disl_A,disl_A_lm,disl_A_lm,disl_A_lm,disl_A_lm,\
               disl_crust_A, disl_A, disl_A_lm, disl_A_lm,\
                disl_A,disl_A,disl_A,disl_A,disl_A_lm,disl_A_lm,disl_A_lm,disl_A_lm,\
                    disl_crust_A, disl_A, disl_A_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Stress exponents for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (disl_n,disl_n,disl_n,disl_n,disl_n_lm,disl_n_lm,disl_n_lm,disl_n_lm,\
                disl_crust_n, disl_n, disl_n_lm, disl_n_lm,\
                      disl_n,disl_n,disl_n,disl_n,disl_n_lm,disl_n_lm,disl_n_lm,disl_n_lm,\
                          disl_crust_n, disl_n, disl_n_lm)
        
        o_dict['Material model']['Visco Plastic TwoD']['Activation energies for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (disl_E,disl_E,disl_E,disl_E,disl_E_lm,disl_E_lm,disl_E_lm,disl_E_lm,\
                disl_crust_E, disl_E, disl_E_lm, disl_E_lm,\
                    disl_E,disl_E,disl_E,disl_E,disl_E_lm,disl_E_lm,disl_E_lm,disl_E_lm,\
                        disl_crust_E, disl_E, disl_E_lm)
    
        o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for dislocation creep'] = \
            "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    sp_upper: %.4e|%.4e|%.4e|%.4e,\
    sp_lower: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
    plate_edge: %.4e|%.4e|%.4e"\
            % (disl_V,disl_V,disl_V,disl_V,disl_V_lm,disl_V_lm,disl_V_lm,disl_V_lm,\
                disl_crust_V, disl_V, disl_V_lm, disl_V_lm,\
                    disl_V,disl_V,disl_V,disl_V,disl_V_lm,disl_V_lm,disl_V_lm,disl_V_lm,\
                        disl_crust_V, disl_V, disl_V_lm)


def wb_configure_transit_ov_plates(i_feature, trench, ov_age,\
    ov_trans_age, ov_trans_length, wb_new_ridge, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
    geometry = kwargs.get('geometry', 'chunk')
    side_angle = 360  # side angle to creat features in the 3rd dimension
    side_dist = 10000e3
    Ro = kwargs.get("Ro", 6371e3)
    o_feature = i_feature.copy()
    trans_angle = ov_trans_length / Ro / np.pi * 180.0
    if geometry == 'chunk':
        ov = trench  + trans_angle  # new ending point of the default overiding plage
        side = side_angle
        ridge = trench - trans_angle * ov_trans_age / (ov_age - ov_trans_age)
    elif geometry == 'box':
        ov = trench + ov_trans_length
        side = side_dist
        ridge = trench - ov_trans_length * ov_trans_age / (ov_age - ov_trans_age)
    else:
        pass
    v = ov_trans_length / (ov_age - ov_trans_age)
    o_feature["temperature models"][0]["spreading velocity"] = v
    o_feature["coordinates"] = [[trench, -side], [trench, side],\
        [ov, side], [ov, -side]]
    if wb_new_ridge == 1:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[[ridge, -side], [ridge, side]]]
    else:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[ridge, -side], [ridge, side]]
    return o_feature, ov


# todo_ov
def wb_configure_transit_ov_plates_1(i_feature, trench, ov_age,\
    ov_trans_age, ov_trans_length, wb_new_ridge, sp_width, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
    geometry = kwargs.get('geometry', 'chunk')
    side_angle = 360  # side angle to creat features in the 3rd dimension
    Ro = kwargs.get("Ro", 6371e3)
    o_feature = i_feature.copy()
    trans_angle = ov_trans_length / Ro / np.pi * 180.0
    if geometry == 'chunk':
        ov = trench  + trans_angle  # new ending point of the default overiding plage
        side = side_angle # I haven't implement for chuck, so I just keep a big value here
        ridge = trench - trans_angle * ov_trans_age / (ov_age - ov_trans_age)
    elif geometry == 'box':
        ov = trench + ov_trans_length
        side = sp_width
        ridge = trench - ov_trans_length * ov_trans_age / (ov_age - ov_trans_age)
    else:
        pass
    v = ov_trans_length / (ov_age - ov_trans_age)
    o_feature["temperature models"][0]["spreading velocity"] = v
    o_feature["coordinates"] = [[trench, -side], [trench, side],\
        [ov, side], [ov, -side]]
    if wb_new_ridge == 1:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[[ridge, -side], [ridge, side]]]
    else:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[ridge, -side], [ridge, side]]
    return o_feature, ov


def wb_configure_transit_ov_plates_1_chunk(i_feature, trench, ov_age,\
    ov_trans_age, ov_trans_length, wb_new_ridge, sp_width, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
    Ro = kwargs.get("Ro", 6371e3)
    trench_deg = trench / Ro * 180.0 / np.pi
    o_feature = i_feature.copy()
    ov = trench  + ov_trans_length  # new ending point of the default overiding plage
    side_deg = sp_width / Ro * 180.0 / np.pi
    ridge = trench - ov_trans_length * ov_trans_age / (ov_age - ov_trans_age)
    ov_deg = ov / Ro * 180.0 / np.pi
    ridge_deg = ridge / Ro * 180.0 / np.pi
    v = ov_trans_length / (ov_age - ov_trans_age)
    o_feature["temperature models"][0]["spreading velocity"] = v
    o_feature["coordinates"] = [[trench_deg, -1.0], [trench_deg, side_deg],\
        [ov_deg, side_deg], [ov_deg, -1.0]]
    if wb_new_ridge == 1:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[[ridge_deg, -1.0], [ridge_deg, side_deg]]]
    else:
        o_feature["temperature models"][0]["ridge coordinates"] =\
            [[ridge_deg, -1.0], [ridge_deg, side_deg]]
    return o_feature, ov


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - create case with json file: \n\
\n\
        Lib_ThDSubduction0_Cases create_with_json -j \
        /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_create_test/configure_1.json"
        )


def ShowJsonOption():
    Case_Opt = CASE_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Case_Opt.document_str()
        )


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='path to a json file')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif (_commend in ['--json_option', '-jo']):
        # json options
        ShowJsonOption()
    elif _commend == 'create_with_json':
        # example:
        CasesP.create_case_with_json(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()