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
    First, substitute FOO with the name of the project
"""
import numpy as np
import sys, os, argparse
# import pathlib
import numpy as np
import shilofue.ParsePrm as ParsePrm
import shilofue.PlotVisit as PlotVisit
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, RunScripts
import shilofue.TwoDSubduction0.PlotVisit as TwoDPlotVisit

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
        Lib_FOO0_PlotVisit visit_options -i `pwd` -sr temperature.py\n\
            -sr: the relative path under the visit_script folder\n\
\n\
  - run script: \n\
        Lib_FOO0_PlotVisit run -i `pwd`\n\
\n\
  - run vtk scripts: \n\
        Lib_FOO0_PlotVisit vtk_options -i `pwd` -p FOO_MOW -s 25\n\
            -p operations, available operations are: \n\
                FOO_MOW: pull out tentative MOW from temperature: \n\
                FOO_SlabAnalysis: analyze slab morphology: \n\
        ")


class VISIT_OPTIONS(TwoDPlotVisit.VISIT_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, **kwargs):
        """
        Interpret the inputs, to be reloaded in children
        kwargs: options
            last_step(list): plot the last few steps
        """
        # additional inputs
        rotation_plus = kwargs.get("rotation_plus", 0.0) # additional rotation
        
        to_rad = np.pi / 180.0
        # call function from parent
        PlotVisit.VISIT_OPTIONS.Interpret(self, **kwargs)
        idx = ParsePrm.FindWBFeatures(self.wb_dict, "Subducting plate")
        idx1 = ParsePrm.FindWBFeatures(self.wb_dict, "Slab")

        # Geometry
        sub_plate_feature = self.wb_dict["features"][idx]
        slab_feature = self.wb_dict["features"][idx1]
        # this is point marking the extension of plate
        sub_plate_extends = sub_plate_feature['coordinates'][1]
        box_width = -1.0
        Ro = -1.0
        sp_age = np.nan
        ov_age = np.nan
        self.options['ROTATION_ANGLE'] = 0.0
        if self.options["GEOMETRY"] == "box":
            box_width = self.idict["Geometry model"]["Box"]["Y extent"]
            self.options['TRENCH_EDGE_Y'] = sub_plate_extends[1] * 0.75
            self.options['TRENCH_EDGE_Y_FULL'] = sub_plate_extends[1]
            self.options['BOX_WIDTH'] = box_width
            self.options['BOX_THICKNESS'] = self.idict["Geometry model"]["Box"]["Z extent"]
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except ParsePrm.WBFeatureNotFoundError:
                pass
            else:
                feature_sp = self.wb_dict['features'][index]
                trench_x = feature_sp["coordinates"][2][0]
                try:
                    spreading_velocity = feature_sp["temperature models"][0]["spreading velocity"]
                except KeyError:
                    pass
                else:
                    sp_age = trench_x / spreading_velocity 
            try: 
                idx2 = ParsePrm.FindWBFeatures(self.wb_dict, "Overiding plate")
            except ParsePrm.WBFeatureNotFoundError:
                pass
            else:
                ov_plate_feature = self.wb_dict["features"][idx2]
                ov_age = ov_plate_feature["temperature models"][0]["plate age"]
            self.options['ROTATION_ANGLE'] = 0.0
        elif self.options["GEOMETRY"] == "chunk":
            # todo_3d_chunk
            # in chunk geometry, the coordinate is read in as the latitude, and it's in
            # degree
            Ro = float(self.idict["Geometry model"]["Chunk"]["Chunk outer radius"])
            self.options['TRENCH_EDGE_Y'] = sub_plate_extends[1] * np.pi / 180.0 * Ro * 0.75
            self.options['TRENCH_EDGE_Y_FULL'] = sub_plate_extends[1] * np.pi / 180.0 * Ro
            self.options['TRENCH_EDGE_LAT_FULL'] = sub_plate_extends[1]
            self.options["CHUNK_RIDGE_CENTER_X"] = Ro
            self.options["CHUNK_RIDGE_CENTER_Z"] = 0.0
            # convert to x, y, z with long = 0.75 * plate_extent, lat = 0.0, r = Ro
            ridge_edge_x, ridge_edge_y, ridge_edge_z = Utilities.ggr2cart(sub_plate_extends[1]*to_rad*0.75, 0.0, Ro)
            self.options["CHUNK_RIDGE_EDGE_X"] = ridge_edge_x
            self.options["CHUNK_RIDGE_EDGE_Z"] = ridge_edge_z
            self.options['BOX_WIDTH'] = -1.0
            try:
                index = ParsePrm.FindWBFeatures(self.wb_dict, 'Subducting plate')
            except ParsePrm.WBFeatureNotFoundError:
                # either there is no wb file found, or the feature 'Subducting plate' is not defined
                rotation_angle = 52.0 + rotation_plus
            else:
                # rotate to center on the slab
                feature_sp = self.wb_dict['features'][index]
                trench_phi = feature_sp["coordinates"][2][0]
                rotation_angle = 90.0 - trench_phi - 2.0 + rotation_plus
                try:
                    spreading_velocity = feature_sp["temperature models"][0]["spreading velocity"]
                except KeyError:
                    pass
                else:
                    sp_age = trench_phi * np.pi / 180.0 * Ro/ spreading_velocity 
            try: 
                idx2 = ParsePrm.FindWBFeatures(self.wb_dict, "Overiding plate")
            except ParsePrm.WBFeatureNotFoundError:
                pass
            else:
                ov_plate_feature = self.wb_dict["features"][idx2]
                ov_age = ov_plate_feature["temperature models"][0]["plate age"]
                self.options['ROTATION_ANGLE'] = rotation_angle
        else:
            raise ValueError("geometry must by either box or chunk")


        # viscosity
        self.options['ETA_MIN'] = self.idict['Material model']['Visco Plastic TwoD']['Minimum viscosity']
        self.options['ETA_MAX'] = self.idict['Material model']['Visco Plastic TwoD']['Maximum viscosity']
        self.options['TRENCH_INITIAL'] = slab_feature['coordinates'][1][0] 
        
        # yield stress
        try:
            self.options["MAXIMUM_YIELD_STRESS"] = float(self.idict['Material model']['Visco Plastic TwoD']["Maximum yield stress"])
        except KeyError:
            self.options["MAXIMUM_YIELD_STRESS"] = 1e9

        # peierls rheology
        try:
            include_peierls_rheology = self.idict['Material model']['Visco Plastic TwoD']['Include Peierls creep']
            if include_peierls_rheology == 'true':
                self.options['INCLUDE_PEIERLS_RHEOLOGY'] = True
            else:
                self.options['INCLUDE_PEIERLS_RHEOLOGY'] = False
        except KeyError:
            self.options['INCLUDE_PEIERLS_RHEOLOGY'] = False

        # age
        self.options["SP_AGE"] = sp_age
        self.options["OV_AGE"] = ov_age

    def vtk_options(self, **kwargs):
        '''
        options of vtk scripts
        '''
        # call function from parent
        PlotVisit.VISIT_OPTIONS.vtk_options(self, **kwargs)
    
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
            center_profile_file_path = os.path.join(self._case_dir, "vtk_outputs", "center_profile_%05d.txt" % snap)
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
