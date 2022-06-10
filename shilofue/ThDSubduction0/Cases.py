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

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

year = 365 * 24 * 3600.0  # yr to s

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

    
    def check(self):
        reset_trailing_morb = self.values[self.start+7]
        assert(reset_trailing_morb in [0, 1])
        friction_angle = self.values[self.start+10] # range of friction angle, in degree
        assert(friction_angle >= 0.0 and friction_angle <= 90.0)

    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
        _type = self.values[9] 
        box_width = self.values[self.start]
        box_length = self.values[self.start+1]
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
        return _type, if_wb, geometry, box_width, box_length, box_depth,\
            sp_width, trailing_length, reset_trailing_morb, ref_visc,\
            relative_visc_plate, friction_angle, relative_visc_lower_mantle, cohesion,\
            sp_depth_refining, reference_density, sp_relative_density, global_refinement,\
            adaptive_refinement
        

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
        return _type, if_wb, geometry, sp_width, sp_length, trailing_length



class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, _type, if_wb, geometry, box_width, box_length, box_depth,\
    sp_width, trailing_length, reset_trailing_morb, ref_visc, relative_visc_plate, friction_angle,\
    relative_visc_lower_mantle, cohesion, sp_depth_refining, reference_density, sp_relative_density, \
    global_refinement, adaptive_refinement):
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
        # rheology
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
        # rewrite the reaction morb part
        if _type == 's07':
            if reset_trailing_morb == 1:
                o_dict['Material model'][material_model_subsection]['Reaction mor'] = 'true'
                o_dict['Material model'][material_model_subsection]['Reaction mor function']['Function constants'] =\
                    "Width=%.4e, Do=%.4e, xm=%.4e, DpUp=5e+04, Dp=1e5, Wp=%.4e, pWidth=1e5" %  (trailing_length, box_depth, box_length, sp_width)
            else:
                o_dict['Material model'][material_model_subsection]['Reaction mor'] = 'false'

        o_dict['Material model'][material_model_subsection] = {**o_dict['Material model'][material_model_subsection], **outputs}  # prepare entries
            

        pass


    def configure_wb(self, _type, if_wb, geometry, sp_width, sp_length, trailing_length):
        '''
        Configure wb file
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # geometry options
        wb_configure_plate_schellart07(self.wb_dict, sp_width, sp_length, trailing_length)

        pass


def wb_configure_plate_schellart07(wb_dict, sp_width, sp_length, trailing_width):
    '''
    World builder configuration of plates in Schellart etal 2007
    '''
    # subducting plate
    o_dict = wb_dict.copy()
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[trailing_width, -sp_width], [trailing_width, sp_width], [sp_length, sp_width] ,[sp_length, -sp_width]]
    o_dict['features'][i0] = sp_dict
    # slab
    o_dict = wb_dict.copy()
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    sdict = o_dict['features'][i0]
    sdict["coordinates"] = [[sp_length, -sp_width], [sp_length, sp_width]]
    o_dict['features'][i0] = sdict





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