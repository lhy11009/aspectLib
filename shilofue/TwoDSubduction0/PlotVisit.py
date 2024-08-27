# -*- coding: utf-8 -*-
r"""Plot visit output

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
# import pathlib
import numpy as np
import shilofue.ParsePrm as ParsePrm
import shilofue.PlotVisit as PlotVisit

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
Visualization with visit\n\
\n\
Examples of usage: \n\
\n\
  - translate script: \n\
\n\
        python -m shilofue.TwoDSubduction0.PlotVisit visit_options -i $TwoDSubduction_DIR/latent_heat_issue/cookbook_latent-heat -sr temperature.py\n\
            -sr: the relative path under the visit_script folder\n\
\n\
  - run script: \n\
        python -m shilofue.TwoDSubduction0.PlotVisit run -i $TwoDSubduction_DIR/non_linear34/eba_low_tol_newton_shift_CFL0.8/visit_scripts/slab.py\n\
\n\
  - run vtk scripts: \n\
        python -m shilofue.TwoDSubduction0.PlotVisit vtk_options -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8_lh -p TwoDSubduction_MOW -s 25\n\
            -p operations, available operations are: \n\
                TwoDSubduction_MOW: pull out tentative MOW from temperature: \n\
                TwoDSubduction_SlabAnalysis: analyze slab morphology: \n\
        ")


class VISIT_OPTIONS(PlotVisit.VISIT_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
            last_step(list): plot the last few steps
        """
        # call function from parent
        PlotVisit.VISIT_OPTIONS.Interpret(self, **kwargs)

        # default settings
        self.options['IF_PLOT_SHALLOW'] = kwargs.get('if_plot_shallow', "False") # plot the shallow part of the slab.
        self.options["PLOT_TYPES"] = str(kwargs.get('plot_types', ["upper_mantle"]))
        self.options['IF_EXPORT_SLAB_MORPH'] = 'False'
        self.options['IF_PLOT_SLAB'] = 'True'
        self.options['GLOBAL_UPPER_MANTLE_VIEW_BOX'] = 0.0
        self.options['ROTATION_ANGLE'] = 0.0

        # additional inputs
        rotation_plus = kwargs.get("rotation_plus", 0.0) # additional rotation

        # try using the value for the background
        try:
            self.options['ETA_MIN'] =\
                Utilities.string2list(self.idict['Material model']['Visco Plastic TwoD']['Minimum viscosity'], float)[0]
        except ValueError:
            eta_min_inputs =\
                ParsePrm.COMPOSITION(self.idict['Material model']['Visco Plastic TwoD']['Minimum viscosity']) 
            self.options['ETA_MIN'] = eta_min_inputs.data['background'][0] # use phases
        try:
            self.options['ETA_MAX'] =\
                Utilities.string2list(self.idict['Material model']['Visco Plastic TwoD']['Maximum viscosity'], float)[0]
        except ValueError:
            eta_max_inputs =\
                ParsePrm.COMPOSITION(self.idict['Material model']['Visco Plastic TwoD']['Maximum viscosity']) 
            self.options['ETA_MAX'] = eta_max_inputs.data['background'][0] # use phases
        # self.options['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
        self.options['IF_DEFORM_MECHANISM'] = 1
        # crustal layers
        # todo_2l
        is_crust_2l = False
        composition_fields = []
        temp_list = self.idict['Compositional fields']['Names of fields'].split(",")
        for temp in temp_list:
            temp1 = Utilities.re_neat_word(temp)
            if temp1 != "":
                composition_fields.append(temp1)
        if "spcrust_up" in composition_fields:
            is_crust_2l = True
        # currently only works for chunk geometry
        if self.options['GEOMETRY'] == 'chunk':
            sp_age = -1.0
            ov_age = -1.0
            Ro = 6371e3
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, "Subducting plate")
                index1 = ParsePrm.FindWBFeatures(self.wb_dict, "Overiding plate")
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                sp_age = -1.0
                ov_age = -1.0
            else:
                feature_sp = self.wb_dict['features'][index]
                feature_ov = self.wb_dict['features'][index1]
                trench_angle = feature_sp["coordinates"][2][0]
                spreading_velocity = feature_sp["temperature models"][0]["spreading velocity"]
                sp_age = trench_angle * np.pi / 180.0 * Ro/ spreading_velocity 
                ov_age = feature_ov["temperature models"][0]["plate age"]
            self.options['SP_AGE'] = sp_age
            self.options['OV_AGE'] =  ov_age
        elif self.options['GEOMETRY'] == 'box':
            sp_age = -1.0
            ov_age = -1.0
            try:
                # todo_ptable
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
                index1 = ParsePrm.FindWBFeatures(self.wb_dict, "Overiding plate")
                feature_sp = self.wb_dict['features'][index]
                feature_ov = self.wb_dict['features'][index1]
                trench_x = feature_sp["coordinates"][2][0]
                spreading_velocity = feature_sp["temperature models"][0]["spreading velocity"]
                sp_age = trench_x / spreading_velocity 
                ov_age = feature_ov["temperature models"][0]["plate age"]
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                sp_age = -1.0
                ov_age = -1.0
                pass
            self.options['SP_AGE'] = sp_age
            self.options['OV_AGE'] =  ov_age
        else:
            raise ValueError("Geometry should be \"chunk\" or \"box\"")
        # rotation of the domain
        if self.options['GEOMETRY'] == 'chunk':
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                rotation_angle = 52.0 + rotation_plus
            else:
                # rotate to center on the slab
                feature_sp = self.wb_dict['features'][index]
                rotation_angle = 90.0 - feature_sp["coordinates"][2][0] - 2.0 + rotation_plus
            self.options['ROTATION_ANGLE'] = rotation_angle
        elif self.options['GEOMETRY'] == 'box':
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                trench_x = 4e6
            else:
                # reset the view point
                feature_sp = self.wb_dict['features'][index]
                trench_x = feature_sp["coordinates"][2][0]
            window_width = 1.8e6
            self.options['GLOBAL_UPPER_MANTLE_VIEW_BOX'] =\
            "(%.4e, %.4e, 1.9e6, 2.9e6)" % (trench_x - window_width/2.0, trench_x + window_width/2.0)
        else:
            raise ValueError("Geometry should be \"chunk\" or \"box\"")
        # Slab configuration
        index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
        feature_sp = self.wb_dict['features'][index]
        # shear zone:
        #   the initial thickness is parsed from the wb file
        #   parse the cutoff depth if the viscosity is decoupled from the eclogite transition
        use_lookup_table_morb = self.idict['Material model']['Visco Plastic TwoD'].get("Use lookup table morb", 'false')
        sz_method = 0
        if use_lookup_table_morb == 'true':
            phase_rheology_mixing_models_str = self.idict['Material model']['Visco Plastic TwoD'].get('Phase rheology mixing models', "0, 0, 0, 0, 0")
            phase_rheology_mixing_models = Utilities.string2list(phase_rheology_mixing_models_str, int)
            sz_method = phase_rheology_mixing_models[1]
        elif use_lookup_table_morb == 'false':
            pass
        else:
            raise ValueError()
        self.options["SHEAR_ZONE_METHOD"] = sz_method
        self.options["INITIAL_SHEAR_ZONE_THICKNESS"] = feature_sp["composition models"][0]["max depth"]
        decoupling_eclogite_viscosity = self.idict['Material model']['Visco Plastic TwoD'].get('Decoupling eclogite viscosity', 'false')
        if decoupling_eclogite_viscosity == 'true':
            self.options["SHEAR_ZONE_CUTOFF_DEPTH"] = float(self.idict['Material model']['Visco Plastic TwoD']["Eclogite decoupled viscosity"]["Decoupled depth"])
        else:
            self.options["SHEAR_ZONE_CUTOFF_DEPTH"] = -1.0
        #  the shear zone constant viscosity is calculated from the prefactor of spcrust
        A_diff_inputs = ParsePrm.COMPOSITION(self.idict['Material model']['Visco Plastic TwoD']['Prefactors for diffusion creep'])
        # todo_2l
        if is_crust_2l:
            self.options["SHEAR_ZONE_CONSTANT_VISCOSITY"] = 1.0 / 2.0 / A_diff_inputs.data['spcrust_up'][0] # use phases
        else:
            self.options["SHEAR_ZONE_CONSTANT_VISCOSITY"] = 1.0 / 2.0 / A_diff_inputs.data['spcrust'][0] # use phases
        # yield stress
        try:
            self.options["MAXIMUM_YIELD_STRESS"] = float(self.idict['Material model']['Visco Plastic TwoD']["Maximum yield stress"])
        except KeyError:
            self.options["MAXIMUM_YIELD_STRESS"] = 1e8

        # peierls rheology
        try:
            include_peierls_rheology = self.idict['Material model']['Visco Plastic TwoD']['Include Peierls creep']
            if include_peierls_rheology == 'true':
                self.options['INCLUDE_PEIERLS_RHEOLOGY'] = True
            else:
                self.options['INCLUDE_PEIERLS_RHEOLOGY'] = False
        except KeyError:
            self.options['INCLUDE_PEIERLS_RHEOLOGY'] = False

    def vtk_options(self, **kwargs):
        '''
        options of vtk scripts
        '''
        # call function from parent
        PlotVisit.VISIT_OPTIONS.vtk_options(self, **kwargs)
        # reference trench point
        self.options['THETA_REF_TRENCH'] = 0.0  # initial value
        if self.options['GEOMETRY'] == 'chunk':
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                theta_ref_trench = 0.63
            else:
                # rotate to center on the slab
                feature_sp = self.wb_dict['features'][index]
                theta_ref_trench = feature_sp["coordinates"][2][0] / 180.0 * np.pi
        elif self.options['GEOMETRY'] == 'box':
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                # for the box geometry, this is the x distance of the trench
                theta_ref_trench = 4000000.0
            else:
                # rotate to center on the slab
                feature_sp = self.wb_dict['features'][index]
                theta_ref_trench = feature_sp["coordinates"][2][0]        
        else:
            raise ValueError("Value of geometry must be either \"chunk\" or \"box\"")
        self.options['THETA_REF_TRENCH'] = theta_ref_trench


    def get_snaps_for_slab_morphology(self, **kwargs):
        '''
        get the snaps for processing slab morphology
        kwargs (dict):
            time_interval (float)
        '''
        ptime_start = kwargs.get('time_start', None)
        ptime_interval = kwargs.get('time_interval', None)
        ptime_end = kwargs.get('time_end', None)
        assert(ptime_interval is None or type(ptime_interval) == float)      
        # steps for processing slab morphology
        snaps, times, _ = PlotVisit.GetSnapsSteps(self._case_dir, 'graphical')
        psnaps = []
        ptimes = []
        last_time = -1e8  # initiate as a small value, so first step is included
        # loop within all the available steps, find steps satisfying the time interval requirement.
        for i in range(len(times)):
            time = times[i]
            snap = snaps[i]
            if ptime_start is not None and time < ptime_start:
                # check on the start
                continue
            if ptime_end is not None and time > ptime_end:
                # check on the end
                break
            if type(ptime_interval) == float:
                if (time - last_time) < ptime_interval:
                    continue  # continue if interval is not reached
            pvtu_file_path = os.path.join(self.options["DATA_OUTPUT_DIR"], "solution", "solution-%05d.pvtu" % snap)
            if os.path.isfile(pvtu_file_path):
                # append if the file is found
                last_time = time
                psnaps.append(snap)
                ptimes.append(time)
        return psnaps

    def get_times_for_slab_morphology(self, **kwargs):
        '''
        get the snaps for processing slab morphology
        kwargs (dict):
            time_interval (float)
        '''
        ptime_start = kwargs.get('time_start', None)
        ptime_interval = kwargs.get('time_interval', None)
        ptime_end = kwargs.get('time_end', None)
        assert(ptime_interval is None or type(ptime_interval) == float)      
        # steps for processing slab morphology
        snaps, times, _ = PlotVisit.GetSnapsSteps(self._case_dir, 'graphical')
        psnaps = []
        ptimes = []
        last_time = -1e8  # initiate as a small value, so first step is included
        # loop within all the available steps, find steps satisfying the time interval requirement.
        for i in range(len(times)):
            time = times[i]
            snap = snaps[i]
            if ptime_start is not None and time < ptime_start:
                # check on the start
                continue
            if ptime_end is not None and time > ptime_end:
                # check on the end
                break
            if type(ptime_interval) == float:
                if (time - last_time) < ptime_interval:
                    continue  # continue if interval is not reached
            pvtu_file_path = os.path.join(self.options["DATA_OUTPUT_DIR"], "solution", "solution-%05d.pvtu" % snap)
            if os.path.isfile(pvtu_file_path):
                # append if the file is found
                last_time = time
                psnaps.append(snap)
                ptimes.append(time)
        return ptimes
    
    def get_snaps_for_slab_morphology_outputs(self, **kwargs):
        '''
        get the snaps for processing slab morphology, look for existing outputs
        kwargs (dict):
            time_interval (float)
        '''
        ptime_start = kwargs.get('time_start', None)
        ptime_interval = kwargs.get('time_interval', None)
        ptime_end = kwargs.get('time_end', None)
        assert(ptime_interval is None or type(ptime_interval) == float)      
        # steps for processing slab morphology
        snaps, times, _ = PlotVisit.GetSnapsSteps(self._case_dir, 'graphical')
        psnaps = []
        ptimes = []
        last_time = -1e8  # initiate as a small value, so first step is included
        # loop within all the available steps, find steps satisfying the time interval requirement.
        for i in range(len(times)):
            time = times[i]
            snap = snaps[i]
            if ptime_start is not None and time < ptime_start:
                # check on the start
                continue
            if ptime_end is not None and time > ptime_end:
                # check on the end
                break
            if type(ptime_interval) == float:
                if (time - last_time) < ptime_interval:
                    continue  # continue if interval is not reached
            center_profile_file_path = os.path.join(self._case_dir, "vtk_outputs", "slab_morph_s%06d.txt" % snap)
            if os.path.isfile(center_profile_file_path):
                # append if the file is found
                last_time = time
                psnaps.append(snap)
                ptimes.append(time)
        return ptimes, psnaps



class PREPARE_RESULT_OPTIONS(PlotVisit.PREPARE_RESULT_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
        """
        # call function from parent
        PlotVisit.PREPARE_RESULT_OPTIONS.Interpret(self, **kwargs)


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
    parser.add_argument('-sr', '--script', type=str,
                        default='temperature.py',
                        help='script')
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-p', '--operation', type=str,
                        default='',
                        help='operation to take')
    parser.add_argument('-s', '--step', type=int,
                        default='0',
                        help='step')
    parser.add_argument('-ls', '--last_step', type=bool,
                        default=False,
                        help='Only plot last step (bool value)')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend in ['-h', '--help']:
        # example:
        Usage()

    elif _commend == 'visit_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        case_dir = arg.inputs
        Visit_Options = VISIT_OPTIONS(case_dir)
        # call function
        Visit_Options.Interpret()
        # ofile = os.path.join('visit_scripts', 'slab_sph.py')
        ofile = os.path.join('visit_scripts', os.path.basename(arg.script))
        visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', arg.script)
        visit_base_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'base.py')  # base.py : base file
        Visit_Options.read_contents(visit_base_script, visit_script)  # this part combines two scripts
        Visit_Options.substitute()  # substitute keys in these combined file with values determined by Interpret() function
        ofile_path = Visit_Options.save(ofile, relative=True)  # save the altered script
        pass
    
    elif _commend == 'vtk_options':
        vtk_option_path, _, _ = PlotVisit.PrepareVTKOptions(VISIT_OPTIONS, arg.inputs, arg.operation, step=arg.step)
        PlotVisit.RunVTKScripts(arg.operation, vtk_option_path)
    
    elif _commend == 'run':
        PlotVisit.RunScripts(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
