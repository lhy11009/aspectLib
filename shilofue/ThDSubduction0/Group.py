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
import shilofue.Group as GroupP
from shilofue.Group import CreateGroup
from shilofue.ThDSubduction0.Cases import CASE, CASE_OPT
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities
import shilofue.ThDSubduction0.VtkPp as ThDVtkPp


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage, create group: \n\
        Lib_ThDSubduction0_Group create_group -j ~/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/TwoDSubduction/test_group/test.json\n"
        )

def ShowJsonOption():
    Group_Opt = GROUP_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Group_Opt.document_str()
        )


class GROUP_OPT(GroupP.GROUP_OPT):
    pass


class GROUP(GroupP.GROUP):
    pass


class CASE_SUMMARY(GroupP.CASE_SUMMARY):
    '''
    Attributes:
        cases: name of cases
        steps: end steps of cases
        times: end times of cases
        wallclocks: running time of cases on the wall clock
        ab_paths: absolution_paths of cases
        t660s: time the slab tip reaches 660 km
    '''
    def __init__(self, **kwargs):
        '''
        initiation
        Inputs:
            kwargs
        '''
        GroupP.CASE_SUMMARY.__init__(self, **kwargs)
        self.geometries = []
        self.attrs.append("geometries")
        self.t660s = []
        self.attrs.append("t660s")
        self.t800s = []
        self.attrs.append("t800s")
        self.t1000s = []
        self.attrs.append("t1000s")
        self.dip660s = []
        self.attrs.append("dip660s")
        self.sp_ages = []
        self.attrs.append("sp_ages")
        self.ov_ages = []
        self.attrs.append("ov_ages")
        self.slab_strs = []
        self.attrs.append("slab_strs")
        self.attrs_to_output += ["geometries", "t660s", "t800s", "t1000s", "dip660s", "sp_ages", "ov_ages", "slab_strs"]

    def Update(self, **kwargs):
        '''
        Update on properties
        Inputs:
            kwargs:
                actions (list): actions to take
        '''
        actions = kwargs.get('actions', [])

        GroupP.CASE_SUMMARY.Update(self, **kwargs)

        if "geometry" in actions:
            # initiate these field
            self.geometries = [np.nan for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_geometry(i)


        if "t660" in actions:
            # initiate these field
            self.t660s = [np.nan for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_t660(i)
        
        if "t800" in actions:
            # initiate these field
            self.t800s = [np.nan for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_t800(i)
        
        if "t1000" in actions:
            # initiate these field
            self.t1000s = [np.nan for i in range(self.n_case)]
            # update on specific properties
            for i in range(self.n_case):
                self.update_t1000(i)

        if "ages" in actions:
            # ov and sp ages
            self.sp_ages = [np.nan for i in range(self.n_case)]
            self.ov_ages = [np.nan for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_plate_ages(i)
        
        if "strength" in actions:
            self.slab_strs = [-1 for i in range(self.n_case)]
            for i in range(self.n_case):
                self.update_slab_strength(i)
        
        # These fields need to be mannualy assigned, so we
        # only initiation a nan value for the first time
        if "dip660" in actions:
            # initiate these field
            if self.dip660s == []:
                self.dip660s = [np.nan for i in range(self.n_case)]
            else:
                pass
    
    def update_t660(self, i):
        '''
        update t660
        '''
        case_dir = self.ab_paths[i]
        time_interval_for_slab_morphology = 1e6
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)

        try: 
            results = ThDVtkPp.CenterProfileAnalyze(case_dir, time_interval_for_slab_morphology)
            t660, step660 = results['t660'], results['step660']
        except FileNotFoundError:
            self.t660s[i] = np.nan
            # raise e
        else:
            if np.isnan(t660):
                self.t660s[i] = np.nan
            else:
                self.t660s[i] = t660 / 1e6

    def update_t800(self, i):
        '''
        update t800
        '''
        case_dir = self.ab_paths[i]
        time_interval_for_slab_morphology = 1e6
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)

        try: 
            results = ThDVtkPp.CenterProfileAnalyze(case_dir, time_interval_for_slab_morphology)
            t800, step800 = results['t800'], results['step800']
        except FileNotFoundError:
            self.t800s[i] = np.nan
            # raise e
        else:
            if np.isnan(t800):
                self.t800s[i] = np.nan
            else:
                self.t800s[i] = t800 / 1e6
    
    def update_t1000(self, i):
        '''
        update t1000
        '''
        case_dir = self.ab_paths[i]
        time_interval_for_slab_morphology = 1e6
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)

        try: 
            results = ThDVtkPp.CenterProfileAnalyze(case_dir, time_interval_for_slab_morphology)
            t1000, step1000 = results['t1000'], results['step1000']
        except FileNotFoundError:
            self.t1000s[i] = np.nan
            # raise e
        else:
            if np.isnan(t1000):
                self.t1000s[i] = np.nan
            else:
                self.t1000s[i] = t1000 / 1e6
                
    def update_plate_ages(self, i):
        '''
        update plate ages
        '''
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.sp_ages[i] = Visit_Options.options["SP_AGE"] / 1e6
            self.ov_ages[i] = Visit_Options.options["OV_AGE"] / 1e6
        except FileNotFoundError:
            self.sp_ages[i] = np.nan
            self.ov_ages[i] = np.nan

    def update_slab_strength(self, i):
        '''
        Update slab strength properties
        ''' 
        case_dir = self.ab_paths[i]
        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.slab_strs[i] = Visit_Options.options["MAXIMUM_YIELD_STRESS"] / 1e6
        except FileNotFoundError:
            self.slab_strs[i] = np.nan

    def write_file_for_paper(self, o_path):
        '''
        Write output file, with a format we like for figure table
        '''
        # append the case name first and the absolute path last
        table_columns = ["sp_ages", "ov_ages", "slab_strs", "t660s", "t800s", "t1000s", "dip660s"]
        attrs_to_output = ["names"] + table_columns
        headers = ["names"] + table_columns

        # write file
        format = o_path.split('.')[-1]
        if format == "txt":
            with open(o_path, 'w') as fout: 
                self.write(fout, attrs_to_output, headers) 
        elif format == "csv":
            with open(o_path, 'w') as fout: 
                self.write_csv(fout, attrs_to_output, headers)
        elif format == "tex":
            with open(o_path, 'w') as fout: 
                self.export_to_latex_table(fout, attrs_to_output, headers)
        else:
            raise ValueError("%s: filename should end with \".txt\" or \".csv\", but the given name is %s" % (Utilities.func_name(), o_path))

        Utilities.my_assert(os.path.isfile(o_path), FileNotFoundError, "%s: file %s is not written successfully" % (Utilities.func_name(), o_path)) 
        print("%s: Write file %s" % (Utilities.func_name(), o_path))



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
                        help='A json file')
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
    elif _commend == 'create_group':
        GroupP.CreateGroup(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()