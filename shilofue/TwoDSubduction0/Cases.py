# -*- coding: utf-8 -*-
r"""class with interfaces for manipulating an aspect case in TwoDSubduction Project

Logistics:
    I mean to define every feature here, under the folder of a project.
    I expect this processes to be different from model to model (e.g. how to change the size of the domain).
    As we don't want to repeat the process of changing prms in the prm file, but to summariz possible settings here.

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
import json, re
# import pathlib
# import subprocess
import warnings
import numpy as np
from copy import deepcopy
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec
import shilofue.Cases as CasesP
import shilofue.ParsePrm as ParsePrm
import shilofue.FlowLaws as flf
from shilofue.Rheology import RHEOLOGY_OPR, ConvertFromAspectInput, STRENGTH_PROFILE, PlotShearZoneStrengh

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
        self.add_key("Age of the subducting plate at trench", float,\
            ['world builder', 'subducting plate','age trench'], 80e6, nick='sp_age_trench')
        self.add_key("Spreading rate of the subducting plate", float,\
            ['world builder', 'subducting plate', 'sp rate'], 0.05, nick='sp_rate')
        self.add_key("Age of the overiding plate", float,\
            ['world builder', "overiding plate", 'age'], 40e6, nick='ov_age')
        self.add_key("Age of the transit overiding plate", float,\
            ['world builder', "overiding plate", "transit", 'age'], -1.0, nick='ov_trans_age')
        self.add_key("Length of the transit overiding plate", float,\
            ['world builder', "overiding plate", "transit", 'length'], 300e3, nick='ov_trans_length')
        self.add_key("Type of boundary condition\n\
            available options in [all free slip, ]", str,\
            ["boundary condition", "model"], "all free slip", nick='type_of_bd')
        self.add_key("Width of the Box", float,\
            ["box width"], 6.783e6, nick='box_width')
        self.add_key("Method to use for prescribing temperature", str,\
         ["prescribe temperature method"], 'function', nick="prescribe_T_method")
        self.add_key("Method to use for adjusting plate age.\n\
        The default method \"by values\" is to assign values of box_width, sp_rate, and sp_age_trench.\n\
        The method \"adjust box width\" is to change box_width\
        by assuming a default box_width for a default sp_age_trench\
         and extend the box for an older sp_age_trench.",\
         str,\
         ["world builder", "plate age method"], 'by values', nick="plate_age_method")
        self.add_key("Include peierls creep", int, ['include peierls creep'], 0, nick='if_peierls')
        self.add_key("Coupling the eclogite phase to shear zone viscosity",\
         int, ["shear zone", 'coupling the eclogite phase to shear zone viscosity'], 0, nick='if_couple_eclogite_viscosity')
        self.add_key("Width of the Box before adjusting for the age of the trench.\
This is used with the option \"adjust box width\" for configuring plate age at the trench.\
This value is the width of box for a default age (i.e. 80Myr), while the width of box for a\
different age will be adjusted.",\
          float, ["world builder", "box width before adjusting"], 6.783e6, nick='box_width_pre_adjust')
        self.add_key("Model to use for mantle phase transitions", str,\
         ["phase transition model"], 'CDPT', nick="phase_model")
        self.add_key("Root directory for lookup tables", str,\
         ["HeFESTo model", "data directory"], '.', nick="HeFESTo_data_dir")
        self.add_key("Cutoff depth for the shear zone rheology",\
          float, ["shear zone", 'cutoff depth'], 100e3, nick='sz_cutoff_depth')
        self.add_key("Adjust the refinement of mesh with the size of the box", int,\
          ["world builder", "adjust mesh with box width"], 0, nick='adjust_mesh_with_width') 
        self.add_key("Thickness of the shear zone / crust", float, ["shear zone", 'thickness'], 7.5e3, nick='Dsz')
        self.add_key("Refinement scheme", str, ["refinement scheme"], "2d", nick='rf_scheme')
        self.add_key("peierls creep scheme", str, ['peierls creep', 'scheme'], "MK10", nick='peierls_scheme')
        self.add_key("peierls creep, create a 2 stage model. I want to do this because including peierls scheme in the\
intiation stage causes the slab to break in the middle",\
         float, ['peierls creep', 'two stage intial time'], -1.0, nick='peierls_two_stage_time')
        self.add_key("mantle rheology", str, ['mantle rheology', 'scheme'], "HK03_wet_mod", nick='mantle_rheology_scheme')
        self.add_key("Scheme for shear zone viscosity", str, ["shear zone", 'viscous scheme'], "constant", nick='sz_viscous_scheme')
        self.add_key("cohesion", float, ['mantle rheology', 'cohesion'], 50e6, nick='cohesion')
        self.add_key("friction", float, ['mantle rheology', 'friction'], 25.0, nick='friction')
        self.add_key("cohesion in the shear zone", float, ['shear zone', 'cohesion'], 10e6, nick='crust_cohesion')
        self.add_key("friction in the shear zone", float, ['shear zone', 'friction'], 2.8624, nick='crust_friction')
        self.add_key("constant viscosity in the shear zone", float, ['shear zone', 'constant viscosity'], 1e20, nick='sz_constant_viscosity')
        self.add_key("use WB new ridge implementation", int, ['world builder', 'use new ridge implementation'], 0, nick='wb_new_ridge')
        

        pass
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        CasesP.CASE_OPT.check(self)
        # geometry options
        Utilities.my_assert(self.values[3] in ['chunk', 'box'], ValueError,\
        "%s: The geometry for TwoDSubduction cases must be \"chunk\" or \"box\"" \
        % Utilities.func_name())
        if self.values[3] == 'box':
            Utilities.my_assert(self.values[8] == 1, ValueError,\
            "%s: When using the box geometry, world builder must be used for initial conditions" \
            % Utilities.func_name())  # use box geometry, wb is mandatory
        # check the setting for adjust box width
        plate_age_method = self.values[self.start + 8] 
        if plate_age_method == 'adjust box width':
            box_width = self.values[self.start + 6]
            if box_width != self.defaults[self.start + 6]:
                warnings.warn("By using \"adjust box width\" method for subduction plate age\
                box width will be automatically adjusted. Thus the given\
                value is not taken.")
            box_width_pre_adjust = self.values[self.start+11]
            sp_age_trench_default = self.defaults[self.start]  # default value for age at trench
            sp_rate_default = self.defaults[self.start+1]  # default value for spreading rate
            Utilities.my_assert(box_width_pre_adjust > sp_age_trench_default * sp_rate_default, ValueError,\
            "For the \"adjust box width\" method to work, the box width before adjusting needs to be wider\
than the multiplication of the default values of \"sp rate\" and \"age trench\"")
        # check the option for refinement
        # todo_affinity
        refinement_level = self.values[15]
        assert(refinement_level in [-1, 9, 10, 11])  # it's either not turned on or one of the options for the total refinement levels
        # check the method to use for phase transition
        phase_model = self.values[self.start + 12]
        Utilities.my_assert( phase_model in ["CDPT", "HeFESTo"], ValueError,\
        "%s: Models to use for phases must by CDPT or HeFESTo" \
        % Utilities.func_name())
        # check the directory for HeFESTo
        o_dir = self.values[2]
        root_level = self.values[7]
        if phase_model == "HeFESTo":  # check the existence of Hefesto files
            HeFESTo_data_dir = self.values[self.start + 13]
            HeFESTo_data_dir_pull_path = os.path.join(o_dir, ".." * (root_level - 1), HeFESTo_data_dir)
            Utilities.my_assert(os.path.isdir(HeFESTo_data_dir_pull_path),\
            FileNotFoundError, "%s is not a directory" % HeFESTo_data_dir_pull_path)
        # assert scheme to use for refinement
        rf_scheme = self.values[self.start + 17]
        assert(rf_scheme in ['2d', '3d coarse'])
        # assert scheme of peierls creep to use
        peierls_scheme = self.values[self.start + 18]
        assert(peierls_scheme in ['MK10', "MK10p"])
        # assert viscous scheme to use
        sz_viscous_scheme = self.values[self.start + 21]
        assert(sz_viscous_scheme in ["stress dependent", "constant"])
        friction = self.values[self.start + 23]
        assert(friction >= 0.0 and friction < 90.0)  # an angle between 0.0 and 90.0
        crust_friction = self.values[self.start + 25]
        assert(crust_friction >= 0.0 and crust_friction < 90.0)  # an angle between 0.0 and 90.0
        sz_constant_viscosity = self.values[self.start + 26]
        assert(sz_constant_viscosity > 0.0)
        wb_new_ridge = self.values[self.start + 27]
        assert(wb_new_ridge in [0, 1])  # use the new ridge implementation or not

    def to_configure_prm(self):
        if_wb = self.values[8]
        type_of_bd = self.values[self.start + 5]
        sp_rate = self.values[self.start + 1]
        ov_age = self.values[self.start + 2]
        potential_T = self.values[4]
        box_width = self.values[self.start + 6]
        geometry = self.values[3]
        prescribe_T_method = self.values[self.start + 7]
        plate_age_method = self.values[self.start + 8] 
        if plate_age_method == 'adjust box width':
            box_width = re_write_geometry_while_assigning_plate_age(
            *self.to_re_write_geometry_pa()
            ) # adjust box width
        if_peierls = self.values[self.start + 9]
        if_couple_eclogite_viscosity = self.values[self.start + 10]
        phase_model = self.values[self.start + 12]
        HeFESTo_data_dir = self.values[self.start + 13]
        root_level = self.values[7]
        HeFESTo_data_dir_relative_path = os.path.join("../"*root_level, HeFESTo_data_dir)
        sz_cutoff_depth = self.values[self.start+14]
        adjust_mesh_with_width = self.values[self.start+15]
        rf_scheme = self.values[self.start + 17]
        peierls_scheme = self.values[self.start + 18]
        peierls_two_stage_time = self.values[self.start + 19]
        mantle_rheology_scheme = self.values[self.start + 20]
        stokes_linear_tolerance = self.values[11]
        end_time = self.values[12]
        refinement_level = self.values[15]
        case_o_dir = self.values[16]
        sz_viscous_scheme = self.values[self.start + 21]
        cohesion = self.values[self.start + 22]
        friction = self.values[self.start + 23]
        crust_cohesion = self.values[self.start + 24]
        crust_friction = self.values[self.start + 25]
        sz_constant_viscosity = self.values[self.start + 26]
        return if_wb, geometry, box_width, type_of_bd, potential_T, sp_rate,\
        ov_age, prescribe_T_method, if_peierls, if_couple_eclogite_viscosity, phase_model,\
        HeFESTo_data_dir_relative_path, sz_cutoff_depth, adjust_mesh_with_width, rf_scheme,\
        peierls_scheme, peierls_two_stage_time, mantle_rheology_scheme, stokes_linear_tolerance, end_time,\
        refinement_level, case_o_dir, sz_viscous_scheme, cohesion, friction, crust_cohesion, crust_friction, sz_constant_viscosity

    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
        potential_T = self.values[4]
        sp_age_trench = self.values[self.start]
        sp_rate = self.values[self.start + 1]
        ov_age = self.values[self.start + 2]
        ov_trans_age = self.values[self.start + 3]
        ov_trans_length = self.values[self.start + 4]
        if self.values[self.start + 3] < 0.0:
            if_ov_trans = False
        else:
            if_ov_trans = True
        is_box_wider = self.is_box_wider()
        Dsz = self.values[self.start + 16]
        wb_new_ridge = self.values[self.start + 27]
        return if_wb, geometry, potential_T, sp_age_trench, sp_rate, ov_age,\
            if_ov_trans, ov_trans_age, ov_trans_length, is_box_wider, Dsz, wb_new_ridge
    
    def to_re_write_geometry_pa(self):
        box_width_pre_adjust = self.values[self.start+11]
        return box_width_pre_adjust, self.defaults[self.start],\
        self.values[self.start], self.values[self.start+1]
    
    def is_box_wider(self):
        '''
        Return whether we should use a box wider than 90 degree in longtitude
        Return:
            True or False
        '''
        box_width = re_write_geometry_while_assigning_plate_age(
            *self.to_re_write_geometry_pa()
            ) # adjust box width
        if box_width > 1e7:
            return True
        else:
            return False
            

class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, if_wb, geometry, box_width, type_of_bd, potential_T,\
    sp_rate, ov_age, prescribe_T_method, if_peierls, if_couple_eclogite_viscosity, phase_model,\
    HeFESTo_data_dir, sz_cutoff_depth, adjust_mesh_with_width, rf_scheme, peierls_scheme,\
    peierls_two_stage_time, mantle_rheology_scheme, stokes_linear_tolerance, end_time,\
    refinement_level, case_o_dir, sz_viscous_scheme, cohesion, friction, crust_cohesion, crust_friction,\
    sz_constant_viscosity):
        Ro = 6371e3
        # velocity boundaries
        if type_of_bd == "all free slip":  # boundary conditions
            if_fs_sides = True  # use free slip on both sides
        else:
            if_fs_sides = False
        o_dict = self.idict.copy()
        # directory to put outputs
        o_dict["Output directory"] = case_o_dir
        # solver schemes
        if abs((stokes_linear_tolerance-0.1)/0.1) > 1e-6:
            # default is negative, thus do nothing
            o_dict["Solver parameters"]["Stokes solver parameters"]["Linear solver tolerance"] = str(stokes_linear_tolerance)
        # time of computation
        if abs((end_time - 60e6)/60e6) > 1e-6:
            o_dict["End time"] = str(end_time)
        # Adiabatic surface temperature
        o_dict["Adiabatic surface temperature"] = str(potential_T)
        # geometry model
        if geometry == 'chunk':
            max_phi = box_width / Ro * 180.0 / np.pi  # extent in term of phi
            o_dict["Geometry model"] = prm_geometry_sph(max_phi, adjust_mesh_with_width=adjust_mesh_with_width)
        elif geometry == 'box':
            o_dict["Geometry model"] = prm_geometry_cart(box_width, adjust_mesh_with_width=adjust_mesh_with_width)
        # refinement
        if refinement_level > 0:
            # these options only take effects when refinement level is positive
            if refinement_level == 9:
                # this is only an option if the input is positive
                # todo_affinity
                o_dict["Mesh refinement"]["Initial global refinement"] = "5"
                o_dict["Mesh refinement"]["Initial adaptive refinement"] = "4"
            elif refinement_level == 10:
                o_dict["Mesh refinement"]["Initial global refinement"] = "5"
                o_dict["Mesh refinement"]["Initial adaptive refinement"] = "5"
                pass
            elif refinement_level == 11:
                o_dict["Mesh refinement"]["Initial global refinement"] = "6"
                o_dict["Mesh refinement"]["Initial adaptive refinement"] = "5"
            o_dict["Mesh refinement"]["Minimum refinement level"] = o_dict["Mesh refinement"]["Initial global refinement"]
            if geometry == 'chunk':
                o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_sph(refinement_level=refinement_level)
            elif geometry == 'box':
                o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_cart(refinement_level=refinement_level)
        else:
            if geometry == 'chunk':
                o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_sph()
            elif geometry == 'box':
                o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_cart()
        # adjust refinement with different schemes, todo_3d_coarse
        if rf_scheme == "3d_coarse":
            pass
        # boundary temperature model
        if geometry == 'chunk':
            o_dict['Boundary temperature model'] = prm_boundary_temperature_sph()
        elif geometry == 'box':
            o_dict['Boundary temperature model'] = prm_boundary_temperature_cart()
        # set up subsection reset viscosity function
        visco_plastic_twoD = self.idict['Material model']['Visco Plastic TwoD']
        if geometry == 'chunk':
            o_dict['Material model']['Visco Plastic TwoD'] =\
              prm_visco_plastic_TwoD_sph(visco_plastic_twoD, max_phi, if_fs_sides=if_fs_sides)
        elif geometry == 'box':
            o_dict['Material model']['Visco Plastic TwoD'] =\
              prm_visco_plastic_TwoD_cart(visco_plastic_twoD, box_width, if_fs_sides=if_fs_sides)
        # set up subsection Prescribed temperatures
        if type_of_bd == "all free slip":
            o_dict["Prescribe internal temperatures"] = "true"
            if geometry == 'chunk':
                if prescribe_T_method == 'plate model':
                    warnings.warn("plate model only works for cartesian model right now, reset to using function")
                o_dict['Prescribed temperatures'] =\
                    prm_prescribed_temperature_sph(max_phi, potential_T, sp_rate, ov_age)
            elif geometry == 'box':
                o_dict["Prescribe internal temperatures"] = "false" # reset this to false as it doesn't work for now
                if prescribe_T_method == 'function':
                    o_dict['Prescribed temperatures'] =\
                        prm_prescribed_temperature_cart(box_width, potential_T, sp_rate, ov_age)
                elif prescribe_T_method == 'plate model':
                    o_dict['Prescribed temperatures'] =\
                        prm_prescribed_temperature_cart_plate_model(box_width, potential_T, sp_rate, ov_age)
        else:
            # remove this feature if otherwise
            pass
        # Material model
        da_file = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', "depth_average.txt")
        assert(os.path.isfile(da_file))
        Operator = RHEOLOGY_OPR()
        # mantle rheology
        Operator.ReadProfile(da_file)
        if mantle_rheology_scheme == "HK03_wet_mod":  # get the type of rheology
            rheology = Operator.MantleRheology_v1(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=20e3,\
    dVdiff=-5.5e-6, dVdisl=0.0, save_profile=1, dAdiff_ratio=0.33333333333, dAdisl_ratio=1.73205080757, save_json=1)
        elif mantle_rheology_scheme == "HK03":
            # in this one, I don't include F because of the issue related to pressure calibration
            rheology = Operator.MantleRheology_v1(rheology=mantle_rheology_scheme, use_effective_strain_rate=False, save_profile=1, save_json=1)
        else:
            # default is to fix F
            rheology = Operator.MantleRheology_v1(rheology=mantle_rheology_scheme, save_profile=1, save_json=1)
        if mantle_rheology_scheme == "HK03_wet_mod" and sz_viscous_scheme == "constant" and\
            abs(sz_constant_viscosity - 1e20)/1e20 < 1e-6:  # assign the rheology
            pass # this is just the default, so skip. Note here we just skip assigning the mantle rheology in the prm
        else:
            CDPT_assign_mantle_rheology(o_dict, rheology, sz_viscous_scheme=sz_viscous_scheme, sz_constant_viscosity=sz_constant_viscosity)
        self.output_files.append(Operator.output_json)
        self.output_files.append(Operator.output_json_aspect)
        self.output_imgs.append(Operator.output_profile) # append plot of initial conition to figures
        # yielding criteria
        if sz_viscous_scheme == "stress dependent":
            CDPT_assign_yielding(o_dict, cohesion, friction, crust_cohesion=crust_cohesion, crust_friction=crust_friction)
        else:
            CDPT_assign_yielding(o_dict, cohesion, friction)
        # append to initial condition output
        if sz_viscous_scheme == "stress dependent":
            plastic_yielding = {}
            plastic_yielding['cohesion'] = crust_cohesion
            plastic_yielding['friction'] = np.tan(crust_friction * np.pi / 180.0)
            plastic_yielding['type'] = 'Coulumb'
            Operator_Sp = STRENGTH_PROFILE()
            rheology_experiment_dislocation = ConvertFromAspectInput(rheology['dislocation_creep'])
            Operator_Sp.SetRheology(disl=rheology_experiment_dislocation, plastic=plastic_yielding)
            fig_path = os.path.join(ASPECT_LAB_DIR, "results", "shear_zone_strength.png")
            PlotShearZoneStrengh(Operator_Sp, fig_path)
            self.output_imgs.append(fig_path)

        # Include peierls rheology
        if if_peierls:
            try:
                temp = o_dict['Material model']['Visco Plastic TwoD']['Peierls fitting parameters']
            except KeyError as e:
                raise KeyError('The options use Peierls rheology by there are missing parameters in the prm file') from e
            o_dict['Material model']['Visco Plastic TwoD']['Include Peierls creep'] = 'true'
            if peierls_scheme in ["MK10", "MK10p"]: 
                Peierls = flf.GetPeierlsApproxVist('MK10')
                o_dict['Material model']['Visco Plastic TwoD']['Peierls glide parameters p'] = str(Peierls['p'])
                o_dict['Material model']['Visco Plastic TwoD']['Peierls glide parameters q'] = str(Peierls['q'])
                o_dict['Material model']['Visco Plastic TwoD']['Stress exponents for Peierls creep'] = str(Peierls['n'])
                o_dict['Material model']['Visco Plastic TwoD']['Peierls stresses'] = '%.4e' % Peierls['sigp0']
                o_dict['Material model']['Visco Plastic TwoD']['Activation energies for Peierls creep'] = '%.4e' % Peierls['E']
                o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for Peierls creep'] = '0.0'
                A = Peierls['A']
                if phase_model == "CDPT":
                    # note that this part contains the different choices of phases
                    # in order to set up for the lower mantle compositions
                    # a future implementation could indicate in the phases which are lower mantle compositions
                    o_dict['Material model']['Visco Plastic TwoD']['Prefactors for Peierls creep'] = \
                    "background: %.4e|%.4e|%.4e|%.4e|1e-31|1e-31|1e-31|1e-31,\
spcrust: %.4e|%.4e|1e-31|1e-31,\
spharz: %.4e|%.4e|%.4e|%.4e|1e-31|1e-31|1e-31|1e-31,\
opcrust: %.4e, opharz: %.4e" % (A, A, A, A, A, A, A, A, A, A, A, A)
                else:
                    pass  # not implemented
                if peierls_scheme == "MK10p":
                    # add p dependence in the peierls stress
                    G0, Gp = flf.GetPeierlsStressPDependence()
                    o_dict['Material model']['Visco Plastic TwoD']['Peierls shear modulus'] = "%.4e" % G0
                    o_dict['Material model']['Visco Plastic TwoD']['Peierls shear modulus derivative'] =  "%.4e" % Gp
                else:
                    # no p dependence, note: sigp = sigp0*(1 + (Gp/G0)*P1), here we want to set Gp = 0
                    o_dict['Material model']['Visco Plastic TwoD']['Peierls shear modulus derivative'] = "0.0"
                    pass
                # todo
        else:
            o_dict['Material model']['Visco Plastic TwoD']['Include Peierls creep'] = 'false'
        # Couple eclogite viscosity
        if if_couple_eclogite_viscosity:
            o_dict['Material model']['Visco Plastic TwoD']["Decoupling eclogite viscosity"] = 'false'
        else:
            o_dict['Material model']['Visco Plastic TwoD']["Decoupling eclogite viscosity"] = 'true'
            o_dict['Material model']['Visco Plastic TwoD']["Eclogite decoupled viscosity"] =\
                {
                    "Decoupled depth": str(sz_cutoff_depth),
                    "Decoupled depth width": '10e3'
                }
        self.idict = o_dict
        # phase model
        if phase_model == "HeFESTo":
            o_dict['Material model']['Visco Plastic TwoD']["Use lookup table"] = 'true'
            o_dict['Material model']['Visco Plastic TwoD']["Lookup table"]["Data directory"] = HeFESTo_data_dir
            pass
        elif phase_model == "CDPT":
            o_dict['Material model']['Visco Plastic TwoD'].pop("Use lookup table", "Foo")
        # create a multi-stage model
        if if_peierls and (peierls_two_stage_time > 0):
            o_dict1 = deepcopy(o_dict)
            # for stage 1
            o_dict['Material model']['Visco Plastic TwoD']['Include Peierls creep'] = 'false'
            o_dict['End time'] = '%.4e' % peierls_two_stage_time
            # for stage 2
            o_dict1['Resume computation'] = 'true'
            self.model_stages = 2
            self.additional_idicts.append(o_dict1)
            pass 

    def configure_wb(self, if_wb, geometry, potential_T, sp_age_trench, sp_rate, ov_ag,\
        if_ov_trans, ov_trans_age, ov_trans_length, is_box_wider, Dsz, wb_new_ridge):
        '''
        Configure world builder file
        Inputs:
            see description of CASE_OPT
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # potential T
        self.wb_dict['potential mantle temperature'] = potential_T
        # geometry
        if geometry == 'chunk':
            self.wb_dict["coordinate system"] = {"model": "spherical", "depth method": "begin segment"}
            if is_box_wider:
                self.wb_dict["cross section"] = [[0, 0], [360.0, 0.0]]
            else:
                self.wb_dict["cross section"] = [[0, 0], [180.0, 0.0]]
        elif geometry == 'box':
            self.wb_dict["coordinate system"] = {"model": "cartesian"}
            if is_box_wider:
                self.wb_dict["cross section"] = [[0, 0], [1e7, 0.0]]
            else:
                self.wb_dict["cross section"] = [[0, 0], [2e7, 0.0]]
        else:
            raise ValueError('%s: geometry must by one of \"chunk\" or \"box\"' % Utilities.func_name())
        # plates
        if geometry == 'chunk':
            if is_box_wider:
                max_sph = 360.0
            else:
                max_sph = 180.0
            Ro = float(self.idict['Geometry model']['Chunk']['Chunk outer radius'])
            # sz_thickness
            self.wb_dict = wb_configure_plates(self.wb_dict, sp_age_trench,\
            sp_rate, ov_ag, wb_new_ridge, Ro=Ro, if_ov_trans=if_ov_trans, ov_trans_age=ov_trans_age,\
            ov_trans_length=ov_trans_length, geometry=geometry, max_sph=max_sph, sz_thickness=Dsz)
        elif geometry == 'box':
            if is_box_wider:
                Xmax = 2e7
            else:
                Xmax = 1e7  # lateral extent of the box
            self.wb_dict = wb_configure_plates(self.wb_dict, sp_age_trench,\
            sp_rate, ov_ag, wb_new_ridge, Xmax=Xmax, if_ov_trans=if_ov_trans, ov_trans_age=ov_trans_age,\
            ov_trans_length=ov_trans_length, geometry=geometry, sz_thickness=Dsz) # plates
        else:
            raise ValueError('%s: geometry must by one of \"chunk\" or \"box\"' % Utilities.func_name())


def wb_configure_plates(wb_dict, sp_age_trench, sp_rate, ov_age, wb_new_ridge, **kwargs):
    '''
    configure plate in world builder
    '''
    Ro = kwargs.get('Ro', 6371e3)
    Xmax = kwargs.get('Xmax', 7e6)
    max_sph = kwargs.get("max_sph", 180.0)
    geometry = kwargs.get('geometry', 'chunk')
    Dsz = kwargs.get("sz_thickness", None)
    o_dict = wb_dict.copy()
    trench_sph = (sp_age_trench * sp_rate / Ro) * 180.0 / np.pi
    trench_cart = sp_age_trench * sp_rate
    max_cart = 2 * Xmax
    side_angle = 5.0  # side angle to creat features in the 3rd dimension
    side_dist = 1e3
    if geometry == 'chunk':
        _side = side_angle
        _max = max_sph
        trench = trench_sph
    elif geometry == 'box':
        _side = side_dist
        _max = max_cart
        trench = trench_cart
        pass
    if wb_new_ridge == 1:
        sp_ridge_coords = [[[0, -_side], [0, _side]]]
    else:
        sp_ridge_coords = [[0, -_side], [0, _side]]
    # Overiding plate
    if_ov_trans = kwargs.get('if_ov_trans', False)  # transit to another age
    if if_ov_trans and ov_age > (1e6 + kwargs['ov_trans_age']):  # only transfer to younger age
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_feature, ov =\
            wb_configure_transit_ov_plates(wb_dict['features'][i0], trench,\
                ov_age, kwargs['ov_trans_age'], kwargs['ov_trans_length'], wb_new_ridge,\
                Ro=Ro, geometry=geometry)
        o_dict['features'][i0] = ov_trans_feature
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except KeyError:
            pass
        o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = trench
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    op_dict = o_dict['features'][i0]
    op_dict["coordinates"] = [[ov, -_side], [ov, _side],\
        [_max, _side], [_max, -_side]] # trench position
    op_dict["temperature models"][0]["plate age"] = ov_age  # age of overiding plate
    o_dict['features'][i0] = op_dict
    # Subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[0.0, -_side], [0.0, _side],\
        [trench, _side], [trench, -_side]] # trench position
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    o_dict['features'][i0] = sp_dict
    # Slab
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    s_dict = o_dict['features'][i0]
    s_dict["coordinates"] = [[trench, -_side], [trench, _side]] 
    s_dict["dip point"] = [_max, 0.0]
    s_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    s_dict["temperature models"][0]["plate velocity"] = sp_rate
    if sp_age_trench > 100e6:
        # in this case, I'll use the plate model
        s_dict["temperature models"][0]["use plate model as reference"] = True
        s_dict["temperature models"][0]["max distance slab top"] = 150e3
    o_dict['features'][i0] = s_dict
    # mantle for substracting adiabat
    i0 = ParsePrm.FindWBFeatures(o_dict, 'mantle to substract')
    m_dict = o_dict['features'][i0]
    m_dict["coordinates"] =[[0.0, -_side], [0.0, _side],\
        [_max, _side], [_max, -_side]]
    o_dict['features'][i0] = m_dict
    return o_dict

def wb_configure_transit_ov_plates(i_feature, trench, ov_age,\
    ov_trans_age, ov_trans_length, wb_new_ridge, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
    geometry = kwargs.get('geometry', 'chunk')
    side_angle = 5.0  # side angle to creat features in the 3rd dimension
    side_dist = 1e3
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


def prm_geometry_sph(max_phi, **kwargs):
    '''
    reset geometry for chunk geometry
    '''
    adjust_mesh_with_width = kwargs.get("adjust_mesh_with_width")
    inner_radius = 3.481e6
    outer_radius = 6.371e6
    if adjust_mesh_with_width:
        longitude_repetitions = int(outer_radius * max_phi / 180.0 * np.pi / (outer_radius - inner_radius))
    else:
        longitude_repetitions = 2
    o_dict = {
        "Model name": "chunk",
        "Chunk": {
            "Chunk inner radius": "3.481e6",\
            "Chunk outer radius": "6.371e6",\
            "Chunk maximum longitude": "%.4e" % max_phi,\
            "Chunk minimum longitude": "0.0",\
            "Longitude repetitions": "%d" % longitude_repetitions
        }
    }
    return o_dict


def prm_minimum_refinement_sph(**kwargs):
    """
    minimum refinement function for spherical geometry
    """
    Ro = kwargs.get('Ro', 6371e3)
    # todo_affinity
    refinement_level = kwargs.get("refinement_level", 10)
    if refinement_level == 9:
        R_UM = 6
        R_LS = 7
        pass
    elif refinement_level == 10:
        R_UM = 6
        R_LS = 8
    elif refinement_level == 11:
        R_UM = 7
        R_LS = 9
    else:
        raise ValueError("Wrong value %d for the \"refinement_level\"" % refinement_level)
    o_dict = {
      "Coordinate system": "spherical",
      "Variable names": "r,phi,t",
      "Function constants": "Ro=%.4e, UM=670e3, DD=100e3" % Ro,
      "Function expression": "((Ro-r<UM)? \\\n                                   ((Ro-r<DD)? %d: %d): 0.0)" % (R_LS, R_UM)
    }
    return o_dict


def prm_minimum_refinement_cart(**kwargs):
    """
    minimum refinement function for cartesian geometry
    """
    Do = kwargs.get('Do', 2890e3)
    # todo_affinity
    refinement_level = kwargs.get("refinement_level", 10)
    if refinement_level == 9:
        R_UM = 6
        R_LS = 7
        pass
    elif refinement_level == 10:
        R_UM = 6
        R_LS = 8
    elif refinement_level == 11:
        R_UM = 7
        R_LS = 9
    else:
        raise ValueError("Wrong value for the \"refinement_level\"")
    o_dict = {
      "Coordinate system": "cartesian",
      "Variable names": "x, y, t",
      "Function constants": "Do=%.4e, UM=670e3, DD=100e3" % Do,
      "Function expression": "((Do-y<UM)? \\\n                                   ((Do-y<DD)? %d: %d): 0.0)" % (R_LS, R_UM)
    }
    return o_dict


def prm_boundary_temperature_sph():
    '''
    boundary temperature model in spherical geometry
    '''
    o_dict = {
        "Fixed temperature boundary indicators": "bottom, top",
        "List of model names": "spherical constant",
        "Spherical constant": {
            "Inner temperature": "3500", 
            "Outer temperature": "273"
        }
    }
    return o_dict


def prm_boundary_temperature_cart():
    '''
    boundary temperature model in cartesian geometry
    '''
    o_dict = {
        "Fixed temperature boundary indicators": "bottom, top",
        "List of model names": "box",
        "Box": {
            "Bottom temperature": "3500",
            "Top temperature": "273"
            }
    }
    return o_dict


def prm_geometry_cart(box_width, **kwargs):
    '''
    reset geometry for box geometry
    '''
    adjust_mesh_with_width = kwargs.get("adjust_mesh_with_width")
    inner_radius = 3.481e6
    outer_radius = 6.371e6
    if adjust_mesh_with_width:
        x_repetitions = int(box_width / (outer_radius - inner_radius))
    else:
        x_repetitions = 2
    o_dict = {
        "Model name": "box",
        "Box": {
            "X extent": "%.4e" % box_width,
            "Y extent": "2.8900e6",
            "X repetitions": "%d" % x_repetitions
        }
    }
    return o_dict


def prm_visco_plastic_TwoD_sph(visco_plastic_twoD, max_phi, **kwargs):
    '''
    reset subsection Visco Plastic TwoD
    Inputs:
        visco_plastic_twoD (dict): inputs for the "subsection Visco Plastic TwoD"
        part in a prm file
        kwargs(dict):
    '''
    o_dict = visco_plastic_twoD.copy()
    if_fs_sides = kwargs.get('if_fs_sides', True)
    if if_fs_sides:
        # use free slip on both sides, set ridges on both sides
        o_dict['Reset viscosity'] = 'true'
        o_dict['Reset viscosity function'] =\
            prm_reset_viscosity_function_sph(max_phi)
        o_dict["Reaction mor"] = 'true'
        o_dict["Reaction mor function"] =\
            prm_reaction_mor_function_sph(max_phi)
    else:
        # remove the related options
        pass
    return o_dict


def prm_visco_plastic_TwoD_cart(visco_plastic_twoD, box_width, **kwargs):
    '''
    reset subsection Visco Plastic TwoD
    Inputs:
        visco_plastic_twoD (dict): inputs for the "subsection Visco Plastic TwoD"
        part in a prm file
        kwargs(dict):
    '''
    o_dict = visco_plastic_twoD.copy()
    if_fs_sides = kwargs.get('if_fs_sides', True)
    if if_fs_sides:
        # use free slip on both sides, set ridges on both sides
        o_dict['Reset viscosity'] = 'true'
        o_dict['Reset viscosity function'] =\
            prm_reset_viscosity_function_cart(box_width)
        o_dict["Reaction mor"] = 'true'
        o_dict["Reaction mor function"] =\
            prm_reaction_mor_function_cart(box_width)
    else:
        # remove the related options
        pass
    return o_dict


def prm_reset_viscosity_function_sph(max_phi):
    '''
    Default setting for Reset viscosity function in spherical geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Coordinate system": "spherical",
        "Variable names": "r, phi",
        "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e, CV=1e20" % max_phi_in_rad,
        "Function expression": "(((r > Ro - Depth) && ((Ro*phi < Width) || (Ro*(PHIM-phi) < Width)))? CV: -1.0)"
      }
    return odict


def prm_reset_viscosity_function_cart(box_width):
    '''
    Default setting for Reset viscosity function in cartesian geometry
    '''
    odict = {
        "Coordinate system": "cartesian",
        "Variable names": "x, y",
        "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e, CV=1e20" % box_width,
        "Function expression": "(((y > Do - Depth) && ((x < Width) || (xm-x < Width)))? CV: -1.0)"
    }
    return odict


def prm_reaction_mor_function_sph(max_phi):
    '''
    Default setting for Reaction mor function in cartesian geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Coordinate system": "spherical",\
        "Variable names": "r, phi",\
        "Function constants": "Width=2.75e5, Ro=6.371e6, PHIM=%.4e, DCS=7.500e+03, DHS=3.520e+04" % max_phi_in_rad,\
        "Function expression": "((r > Ro - DCS) && (Ro*phi < Width)) ? 0:\\\n                                        ((r < Ro - DCS) && (r > Ro - DHS) && (Ro*phi < Width)) ? 1:\\\n                                        ((r > Ro - DCS) && (Ro*(PHIM - phi) < Width)) ? 2:\\\n                                        ((r < Ro - DCS) && (r > Ro - DHS) && (Ro*(PHIM - phi) < Width)) ? 3: -1"
      }
    return odict


def prm_reaction_mor_function_cart(box_width):
    '''
    Default setting for Reaction mor function in cartesian geometry
    '''
    odict = {
        "Coordinate system": "cartesian",
        "Variable names": "x, y",
        "Function constants": "Width=2.75e5, Do=2.890e6, xm=%.4e, DCS=7.500e+03, DHS=3.520e+04" % box_width,
        "Function expression": "((y > Do - DCS) && (x < Width)) ? 0:\\\n                                        ((y < Do - DCS) && (y > Do - DHS) && (x < Width)) ? 1:\\\n                                        ((y > Do - DCS) && (xm - x < Width)) ? 2:\\\n                                        ((y < Do - DCS) && (y > Do - DHS) && (xm - x < Width)) ? 3: -1"
    }
    return odict


def prm_prescribed_temperature_sph(max_phi, potential_T, sp_rate, ov_age):
    '''
    Default setting for Prescribed temperatures in spherical geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Indicator function": {
          "Coordinate system": "spherical",\
          "Variable names": "r, phi",\
          "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e" % max_phi_in_rad,\
          "Function expression": "(((r>Ro-Depth)&&((r*phi<Width)||(r*(PHIM-phi)<Width))) ? 1:0)"\
        },\
        "Temperature function": {
          "Coordinate system": "spherical",\
          "Variable names": "r, phi",\
          "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e,\\\n                             AGEOP=%.4e, TS=2.730e+02, TM=%.4e, K=1.000e-06, VSUB=%.4e, PHILIM=1e-6" %\
              (max_phi_in_rad, ov_age * year, potential_T, sp_rate / year),\
          "Function expression": "((r*(PHIM-phi)<Width) ? TS+(TM-TS)*(1-erfc(abs(Ro-r)/(2*sqrt(K*AGEOP)))):\\\n\t(phi > PHILIM)? (TS+(TM-TS)*(1-erfc(abs(Ro-r)/(2*sqrt((K*Ro*phi)/VSUB))))): TM)"
        }
    }
    return odict


def prm_prescribed_temperature_cart(box_width, potential_T, sp_rate, ov_age):
    '''
    Default setting for Prescribed temperatures in cartesian geometry
    '''
    odict = {
        "Indicator function": {
          "Coordinate system": "cartesian",
          "Variable names": "x, y",
          "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e" % box_width,
          "Function expression": "(((y>Do-Depth)&&((x<Width)||(xm-x<Width))) ? 1:0)"
        },
        "Temperature function": {
          "Coordinate system": "cartesian",
          "Variable names": "x, y",
          "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e,\\\n                             AGEOP=%.4e, TS=2.730e+02, TM=%.4e, K=1.000e-06, VSUB=%.4e, XLIM=6" %\
                (box_width, ov_age * year, potential_T, sp_rate / year),\
          "Function expression": "(xm-x<Width) ? TS+(TM-TS)*(1-erfc(abs(Do-y)/(2*sqrt(K*AGEOP)))):\\\n\t((x > XLIM)? (TS+(TM-TS)*(1-erfc(abs(Do-y)/(2*sqrt((K*x)/VSUB))))): TM)"
        }
    }
    return odict


def prm_prescribed_temperature_cart_plate_model(box_width, potential_T, sp_rate, ov_age):
    '''
    Default setting for Prescribed temperatures in cartesian geometry using the plate model
    '''
    odict = {
        "Model name": "plate model",
        "Indicator function": {
          "Coordinate system": "cartesian",
          "Variable names": "x, y",
          "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e" % box_width,
          "Function expression": "(((y>Do-Depth)&&((x<Width)||(xm-x<Width))) ? 1:0)"
        },
        "Plate model":{
            "Subducting plate velocity" : "%.4e" % (sp_rate/year)
        }
    }
    return odict


def re_write_geometry_while_assigning_plate_age(box_width0, sp_age0, sp_age, sp_rate):
    '''
    adjust box width with assigned spreading rate of subducting plate and subducting plate age
    Inputs:
        box_width0: default box width
        sp_age0: default plate age
        sp_age: plate age
        sp_rate: spreading rate of the subducting plate
    '''
    box_width = box_width0 + (sp_age - sp_age0) * sp_rate
    return box_width


def CDPT_assign_mantle_rheology(o_dict, rheology, **kwargs):
    '''
    Assign mantle rheology in the CDPT model
    ''' 
    diffusion_creep = rheology['diffusion_creep']
    dislocation_creep = rheology['dislocation_creep']
    diffusion_creep_lm = rheology['diffusion_lm']
    sz_viscous_scheme = kwargs.get("sz_viscous_scheme", "constant")
    sz_constant_viscosity = kwargs.get("sz_constant_viscosity", 1e20)
    if sz_viscous_scheme == "constant":
        diff_crust_A = 1.0 / 2.0 / sz_constant_viscosity
        diff_crust_m = 0.0
        diff_crust_E = 0.0
        diff_crust_V = 0.0
        disl_crust_A = 5e-32
        disl_crust_n = 1.0
        disl_crust_E = 0.0
        disl_crust_V = 0.0
    elif sz_viscous_scheme == "stress dependent":
        diff_crust_A = 5e-32
        diff_crust_m = 0.0
        diff_crust_E = 0.0
        diff_crust_V = 0.0
        disl_crust_A = dislocation_creep['A']
        disl_crust_n = dislocation_creep['n']
        disl_crust_E = dislocation_creep['E']
        disl_crust_V = dislocation_creep['V']
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
    o_dict['Material model']['Visco Plastic TwoD']['Prefactors for diffusion creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
spcrust: %.4e|%.4e|%.4e|%.4e,\
spharz: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
opcrust: %.4e, opharz: %.4e" % (diff_A, diff_A, diff_A, diff_A, diff_A_lm, diff_A_lm,\
diff_A_lm, diff_A_lm, diff_crust_A, diff_A, diff_A_lm, diff_A_lm,\
diff_A, diff_A, diff_A, diff_A, diff_A_lm, diff_A_lm, diff_A_lm, diff_A_lm,\
diff_A, diff_A)
    o_dict['Material model']['Visco Plastic TwoD']['Grain size exponents for diffusion creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
spcrust: %.4e|%.4e|%.4e|%.4e,\
spharz: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
opcrust: %.4e, opharz: %.4e" % (diff_m, diff_m, diff_m, diff_m, diff_m_lm, diff_m_lm,\
diff_m_lm, diff_m_lm, diff_crust_m, diff_m, diff_m_lm, diff_m_lm,\
diff_m, diff_m, diff_m, diff_m, diff_m_lm, diff_m_lm, diff_m_lm, diff_m_lm,\
diff_m, diff_m)
    o_dict['Material model']['Visco Plastic TwoD']['Activation energies for diffusion creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
spcrust: %.4e|%.4e|%.4e|%.4e,\
spharz: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
opcrust: %.4e, opharz: %.4e" % (diff_E, diff_E, diff_E, diff_E, diff_E_lm, diff_E_lm,\
diff_E_lm, diff_E_lm, diff_crust_E, diff_E, diff_E_lm, diff_E_lm,\
diff_E, diff_E, diff_E, diff_E, diff_E_lm, diff_E_lm, diff_E_lm, diff_E_lm,\
diff_E, diff_E)
    o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for diffusion creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
spcrust: %.4e|%.4e|%.4e|%.4e,\
spharz: %.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e,\
opcrust: %.4e, opharz: %.4e" % (diff_V, diff_V, diff_V, diff_V, diff_V_lm, diff_V_lm,\
diff_V_lm, diff_V_lm, diff_crust_V, diff_V, diff_V_lm, diff_V_lm,\
diff_V, diff_V, diff_V, diff_V, diff_V_lm, diff_V_lm, diff_V_lm, diff_V_lm,\
diff_V, diff_V)
    o_dict['Material model']['Visco Plastic TwoD']['Prefactors for dislocation creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|5.0000e-32|5.0000e-32|5.0000e-32|5.0000e-32,\
spcrust: %.4e|%.4e|5.0000e-32|5.0000e-32,\
spharz: %.4e|%.4e|%.4e|%.4e|5.0000e-32|5.0000e-32|5.0000e-32|5.0000e-32,\
opcrust: %.4e, opharz: %.4e" % (disl_A, disl_A, disl_A, disl_A,\
disl_crust_A, disl_A, disl_A, disl_A, disl_A, disl_A, disl_A, disl_A)
    o_dict['Material model']['Visco Plastic TwoD']['Stress exponents for dislocation creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00,\
spcrust: %.4e|%.4e|1.0000e+00|1.0000e+00,\
spharz: %.4e|%.4e|%.4e|%.4e|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00,\
opcrust: %.4e, opharz: %.4e" % (disl_n, disl_n, disl_n, disl_n,\
disl_crust_n, disl_n, disl_n, disl_n, disl_n, disl_n, disl_n, disl_n)
    o_dict['Material model']['Visco Plastic TwoD']['Activation energies for dislocation creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|0.0000e+00|0.0000e+00|0.0000e+00|0.0000e+00,\
spcrust: %.4e|%.4e|0.0000e+00|0.0000e+00,\
spharz: %.4e|%.4e|%.4e|%.4e|0.0000e+00|0.0000e+00|0.0000e+00|0.0000e+00,\
opcrust: %.4e, opharz: %.4e" % (disl_E, disl_E, disl_E, disl_E,\
disl_crust_E, disl_E, disl_E, disl_E, disl_E, disl_E, disl_E, disl_E)
    o_dict['Material model']['Visco Plastic TwoD']['Activation volumes for dislocation creep'] = \
        "background: %.4e|%.4e|%.4e|%.4e|0.0000e+00|0.0000e+00|0.0000e+00|0.0000e+00,\
spcrust: %.4e|%.4e|0.0000e+00|0.0000e+00,\
spharz: %.4e|%.4e|%.4e|%.4e|0.0000e+00|0.0000e+00|0.0000e+00|0.0000e+00,\
opcrust: %.4e, opharz: %.4e" % (disl_V, disl_V, disl_V, disl_V,\
disl_crust_V, disl_V, disl_V, disl_V, disl_V, disl_V, disl_V, disl_V)


def CDPT_assign_yielding(o_dict, cohesion, friction, **kwargs):
    '''
    Assign mantle rheology in the CDPT model
    ''' 
    crust_cohesion = kwargs.get("crust_cohesion", cohesion)
    crust_friction = kwargs.get("crust_friction", friction)
    if abs(cohesion  - 50e6)/50e6 < 1e-6 and abs(friction - 25.0)/25.0 < 1e-6\
    and abs(crust_cohesion  - 50e6)/50e6 < 1e-6 and  abs(crust_friction - 25.0)/25.0 < 1e-6:
        pass  # default conditions
    else:
        o_dict['Material model']['Visco Plastic TwoD']["Angles of internal friction"] = "background: %.4e, spcrust: %.4e, spharz: %.4e, opcrust: %.4e, opharz: %.4e" \
        % (friction, crust_friction, friction, friction, friction)
        o_dict['Material model']['Visco Plastic TwoD']["Cohesions"] = "background: %.4e, spcrust: %.4e, spharz: %.4e, opcrust: %.4e, opharz: %.4e" \
        % (cohesion, crust_cohesion, cohesion, cohesion, cohesion)




def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - create case with json file: \n\
\n\
        Lib_TwoDSubduction0_Cases create_with_json -j \
        /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_create_test/configure_1.json"
        )


def ShowJsonOption():
    Case_Opt = CASE_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Case_Opt.document_str()
        )


def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
    pass


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