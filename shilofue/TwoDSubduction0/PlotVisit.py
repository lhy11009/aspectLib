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
import shilofue.PlotVisit as PlotVisit, PrepareVTKOptions, RunVTKScripts, RunScripts

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
        Lib_PlotVisit visit_options -i $TwoDSubduction_DIR/latent_heat_issue/cookbook_latent-heat -sr temperature.py\n\
            -sr: the relative path under the visit_script folder\n\
\n\
  - run script: \n\
        Lib_PlotVisit run -i $TwoDSubduction_DIR/non_linear34/eba_low_tol_newton_shift_CFL0.8/visit_scripts/slab.py\n\
\n\
  - run vtk scripts: \n\
        Lib_PlotVisit vtk_options -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8_lh -p TwoDSubduction_MOW -s 25\n\
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
        self.options['IF_EXPORT_SLAB_MORPH'] = 'False'
        self.options['IF_PLOT_SLAB'] = 'True'
        self.options['ETA_MIN'] = self.idict['Material model']['Visco Plastic TwoD']['Minimum viscosity']
        self.options['GLOBAL_UPPER_MANTLE_VIEW_BOX'] = 0.0
        self.options['ROTATION_ANGLE'] = 0.0
        # try using the value for the background
        try:
            self.options['ETA_MAX'] =\
                Utilities.string2list(self.idict['Material model']['Visco Plastic TwoD']['Maximum viscosity'], float)[0]
        except ValueError:
            eta_max_inputs =\
                ParsePrm.COMPOSITION(self.idict['Material model']['Visco Plastic TwoD']['Maximum viscosity']) 
            self.options['ETA_MAX'] = eta_max_inputs.data['background'][0] # use phases
        # self.options['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
        self.options['IF_DEFORM_MECHANISM'] = 1
        # rotation of the domain
        if self.options['GEOMETRY'] == 'chunk':
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except KeyError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                rotation_angle = 52.0
            else:
                # rotate to center on the slab
                feature_sp = self.wb_dict['features'][index]
                rotation_angle = 90.0 - feature_sp["coordinates"][2][0] - 2.0
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
        self.options['THETA_REF_TRENCH'] = theta_ref_trench


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
        vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, arg.inputs, arg.operation, step=arg.step)
        RunVTKScripts(arg.operation, vtk_option_path)
    
    elif _commend == 'run':
        RunScripts(arg.inputs)

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()