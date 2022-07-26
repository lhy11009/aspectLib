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
import shilofue.Cases as CasesP
import shilofue.ParsePrm as ParsePrm
from shilofue.Rheology import RHEOLOGY_OPR
from shilofue.TwoDSubduction0.Cases import re_write_geometry_while_assigning_plate_age
import json

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
        self.add_key("Depth of the box", float, ['geometry setup', 'box depth'], 1000e3, nick='box_depth')
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
        self.add_key("mantle rheology", str, ['mantle rheology', 'scheme'], "HK03_wet_mod", nick='mantle_rheology_scheme')
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

    
    def check(self):
        reset_trailing_morb = self.values[self.start+7]
        assert(reset_trailing_morb in [0, 1])
        friction_angle = self.values[self.start+10] # range of friction angle, in degree
        assert(friction_angle >= 0.0 and friction_angle <= 90.0)
        dip_angle = self.values[self.start+23] # initial dipping angle of the slab
        assert(dip_angle > 0 and dip_angle <= 90.0) # an angle between 0 and 90
        setup_method = self.values[self.start+26] # method of seting up slabs
        assert(setup_method in ['manual', '2d_consistent'])


    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
        _type = self.values[9] 
        setup_method = self.values[self.start+26] # method of seting up slabs
        box_width = self.values[self.start]
        box_length = self.values[self.start+1]
        if setup_method == '2d_consistent':
            # adjust box width with the age and plate spreading rate
            box_length = re_write_geometry_while_assigning_plate_age(
            *self.to_re_write_geometry_pa()
            ) # adjust box width
        box_depth = self.values[self.start+2]
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
        return _type, if_wb, geometry, box_width, box_length, box_depth,\
            sp_width, trailing_length, reset_trailing_morb, ref_visc,\
            relative_visc_plate, friction_angle, relative_visc_lower_mantle, cohesion,\
            sp_depth_refining, reference_density, sp_relative_density, global_refinement,\
            adaptive_refinement, mantle_rheology_scheme, Dsz, apply_reference_density, Ddl, reset_trailing_ov_viscosity
        
    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
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
        return _type, if_wb, geometry, sp_width, sp_length, trailing_length, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age, setup_method, sp_rate
    
    def to_re_write_geometry_pa(self):
        '''
        Interface to re_write_geometry_while_assigning_plate_age
        '''
        box_length_pre_adjust = self.values[self.start+28]
        default_sp_age_trench = self.defaults[self.start+24]
        sp_age_trench = self.values[self.start+24]
        sp_rate = self.values[self.start+27] # method of seting up slabs
        return box_length_pre_adjust, default_sp_age_trench, sp_age_trench, sp_rate


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, _type, if_wb, geometry, box_width, box_length, box_depth,\
    sp_width, trailing_length, reset_trailing_morb, ref_visc, relative_visc_plate, friction_angle,\
    relative_visc_lower_mantle, cohesion, sp_depth_refining, reference_density, sp_relative_density, \
    global_refinement, adaptive_refinement, mantle_rheology_scheme, Dsz, apply_reference_density, Ddl,\
    reset_trailing_ov_viscosity):
        '''
        Configure prm file
        '''
        o_dict = self.idict.copy()
        # geometry options
        # repitition, figure this out by deviding the dimensions with a unit value of repitition_slice
        repitition_slice = np.min(np.array([box_length, box_width, box_depth]))  # take the min as the repitition_slice
        o_dict['Geometry model']['Box']['X extent'] =  str(box_length)
        o_dict['Geometry model']['Box']['X repetitions'] = str(int(box_length//repitition_slice))
        o_dict['Geometry model']['Box']['Y extent'] =  str(box_width)
        o_dict['Geometry model']['Box']['Y repetitions'] = str(int(box_width//repitition_slice))
        o_dict['Geometry model']['Box']['Z extent'] =  str(box_depth)
        o_dict['Geometry model']['Box']['Z repetitions'] = str(int(box_depth//repitition_slice))
        
        # refinement
        o_dict["Mesh refinement"]["Initial global refinement"] = str(global_refinement)
        o_dict["Mesh refinement"]["Minimum refinement level"] = str(global_refinement)
        o_dict["Mesh refinement"]["Initial adaptive refinement"] = str(adaptive_refinement)
        # Minimum refinement function
        max_refinement = global_refinement + adaptive_refinement + 1
        if _type == "s07":
            if (abs(sp_depth_refining - 200e3)/200e3 > 1e-6):
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3" % (box_depth, sp_depth_refining)
            else:
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=200e3, Dp=100e3" % (box_depth)
        elif _type == "s07T":
                o_dict["Mesh refinement"]["Minimum refinement function"]["Function constants"] =\
                "Do=%.4e, UM=670e3, DD=%.4e, Dp=100e3, Rd=%d, Rum=%d" % (box_depth, sp_depth_refining, max_refinement-1, max_refinement-2)

        
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
            da_file = os.path.join(ASPECT_LAB_DIR, 'files', 'ThDSubduction', "depth_average.txt")
            assert(os.path.isfile(da_file))
            Operator = RHEOLOGY_OPR()
            # read profile
            Operator.ReadProfile(da_file)
            if mantle_rheology_scheme == "HK03":
                use_effective_strain_rate = False
            else:
                use_effective_strain_rate = True
            rheology = Operator.MantleRheology_v0(rheology=mantle_rheology_scheme,\
            save_profile=1, save_json=1, use_effective_strain_rate=use_effective_strain_rate)
            s07T_assign_mantle_rheology(o_dict, rheology)    
            self.output_files.append(Operator.output_json)
            self.output_files.append(Operator.output_json_aspect)
            self.output_imgs.append(Operator.output_profile) # append plot of initial conition to figures
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
        # rewrite the reset viscosity part
        if _type == 's07':
            o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
            "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CV=1e20, Wp=%.4e" % (trailing_length, box_depth, box_length, sp_width)
        if _type == 's07_newton' and reset_trailing_ov_viscosity == 1:
            o_dict['Material model'][material_model_subsection]['Reset viscosity function']['Function constants'] =\
            "Depth=1.45e5, Width=%.4e, Do=%.4e, xm=%.4e, CV=1e20, Wp=%.4e" % (trailing_length, box_depth, box_length, sp_width)
        # rewrite the reaction morb part
        if _type in ['s07', 's07_newton']:
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
                if _type == "s07":
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                        "Width=%.4e, Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5" %  (trailing_length, box_depth, box_length, Dsz_str, Dp_str, sp_width)
                elif _type == "s07_newton":
                    o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                        "Do=%.4e, xm=%.4e, DpUp=%s, Dp=%s, Wp=%.4e, pWidth=1e5" %  (box_depth, box_length, Dsz_str, Dp_str, sp_width)
                else:
                    raise ValueError()
            else:
                o_dict['Material model'][material_model_subsection]['Reaction mor'] = 'false'

        o_dict['Material model'][material_model_subsection] = {**o_dict['Material model'][material_model_subsection], **outputs}  # prepare entries
        pass


    def configure_wb(self, _type, if_wb, geometry, sp_width, sp_length, trailing_length, Dsz, Ddl, slab_length,\
    dip_angle, sp_age_trench, ov_age, setup_method, sp_rate):
        '''
        Configure wb file
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # geometry options
        if setup_method == 'manual':
            if _type in ["s07", "s07T"]:
                wb_configure_plate_schellart07(self.wb_dict, sp_width, sp_length, trailing_length, Dsz, Ddl, slab_length)
            elif _type in ["s07_newton"]:
                wb_configure_plate_schellart07_Tdependent(self.wb_dict, sp_width, sp_length, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age)
            else:
                raise ValueError("Wrong value for \"type\"")
        elif setup_method == '2d_consistent':
            self.wb_dict = wb_configure_plate_2d_consistent(self.wb_dict, sp_width, sp_rate, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age)
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


def wb_configure_plate_2d_consistent(wb_dict, sp_width, sp_rate, Dsz, Ddl, slab_length, dip_angle, sp_age_trench, ov_age):
    '''
    World builder configuration of plates in Schellart etal 2007
    '''
    Xmax = 40000e3
    pe_width = 55e3 # width of the plate edges
    sp_length = sp_rate * sp_age_trench
    o_dict = wb_dict.copy()
    # cross section
    o_dict["cross section"] = [[0.0, 0.0], [Xmax, 0.0]]
    # overiding plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    ov_dict = o_dict['features'][i0]
    ov_dict["coordinates"] = [[sp_length, -sp_width], [sp_length, sp_width], [Xmax, sp_width] ,[Xmax, -sp_width]]
    ov_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ov_dict
    # overiding plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate edge')
    ove_dict = o_dict['features'][i0]
    ove_dict["coordinates"] = [[sp_length, sp_width], [sp_length, sp_width + pe_width], [Xmax, sp_width + pe_width] ,[Xmax, sp_width]]
    ove_dict["temperature models"][0]["plate age"] = ov_age
    o_dict['features'][i0] = ove_dict
    # subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[0.0, -sp_width], [0.0, sp_width], [sp_length, sp_width] ,[sp_length, -sp_width]]
    sp_dict["composition models"][0]["max depth"] = Dsz  # moho
    sp_dict["composition models"][1]["min depth"] = Dsz
    sp_dict["composition models"][1]["max depth"] = Dsz + Ddl  # lithosphere - athenosphere boundary
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    o_dict['features'][i0] = sp_dict
    # subducting plate edge
    i0 = ParsePrm.FindWBFeatures(o_dict, "Subducting plate edge")
    spe_dict = o_dict['features'][i0]
    spe_dict["coordinates"] = [[0.0, sp_width], [0.0, sp_width + pe_width], [sp_length, sp_width + pe_width] ,[sp_length, sp_width]]
    spe_dict["temperature models"][0]["spreading velocity"] = sp_rate
    o_dict['features'][i0] = spe_dict
    # slab
    with open(twod_default_wb_file, 'r') as fin:
        twod_default_dict = json.load(fin)
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    i1 = ParsePrm.FindWBFeatures(twod_default_dict, 'Slab')
    o_dict['features'][i0] = twod_default_dict['features'][i1].copy() # default options
    sdict = o_dict['features'][i0]  # modify the properties
    sdict["coordinates"] = [[sp_length, -sp_width], [sp_length, sp_width]]
    sdict["dip point"] = [Xmax, 0.0]
    for i in range(len(sdict["segments"])-1):
        # slab compostion, the last one is a phantom for temperature tapering
        segment = sdict["segments"][i]
        segment["composition models"][0]["max distance slab top"] = Dsz
        segment["composition models"][1]["min distance slab top"] = Dsz
        segment["composition models"][1]["max distance slab top"] = Dsz + Ddl
        sdict["segments"][i] = segment
    sdict["temperature models"][0]["ridge coordinates"] = \
        [[0,-10000000.0], [0, 10000000.0]]
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