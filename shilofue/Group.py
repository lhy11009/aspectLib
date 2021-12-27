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
import json, re
from shilofue.Cases import create_case_with_json
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shutil import rmtree, copy

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \
        ")


class FEATURE_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work a feature with GROUP
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Name of the feature", str, ["name"], "foo", nick='name')
        self.add_key("Key in the json file", list, ["key"], [], nick='key')
        self.add_key("Values to apply", list, ["values"], [], nick='values')
        self.add_key("Abbrevation in case name", list, ["abbreviating value options"], [None, None], nick='abbrev_value_options')
        self.add_key("If abbrevate case name by value defined in this featue", int, ["abbreviation by value"], 0, nick='is_abbrev_value')
        self.add_key("Abbrevation by a given string", list, ["abbreviating strings"], [], nick='abbrev_strings')

    def check(self):
        assert(len(self.values[3]) == 2)  # abbreviation has 2 entries (name, scale)
        if self.if_abbrev_value():
            # contains a prefix and a scaling
            assert(len(self.values[3]) == 2)
        else:
            # number of strings equal number of values
            assert(len(self.values[5]) == len(self.values[2]))
        pass

    def get_keys(self):
        return self.values[1]

    def get_values(self):
        return self.values[2]

    def if_abbrev_value(self):
        return self.values[4]
    
    def get_abbrev_value_options(self):
        return self.values[3]

    def get_abbrev_strings(self):
        return self.values[5]


class GROUP_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with GROUP
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Base name", str, ["base name"], "foo", nick='bname')
        self.add_features('Features to look for', ['features'], FEATURE_OPT)
        self.add_key("Base json file", str, ["base json"], "foo.json", nick='base_json')
        self.add_key("Directory to output to", str, ["output directory"], ".", nick='output_dir')
        self.add_key("Bindings in feature", list, ["bindings"], [], nick='bindings')
        self.add_key("Base directory to import", str, ["base directory"], ".", nick='base_dir')
    
    def check(self):
        if self.values[4] != []:
            for binding in self.values[4]:
                assert(len(binding) == len(self.values[1]))  # binding has length of features
                for i in range(len(binding)):
                    assert(type(binding[i]) == int)

    def get_base_json_path(self):
        return Utilities.var_subs(self.values[2])
    
    def get_features(self):
        '''
        return all the contained features
        '''
        return self.values[1]
    
    def get_output_dir(self):
        '''
        return output directory
        '''
        return Utilities.var_subs(self.values[3])
    
    def to_create_group(self):
        return self.values[1], Utilities.var_subs(self.values[5]), Utilities.var_subs(self.values[3]), self.values[0], self.values[4]


class GROUP():
    '''
    Define a class for groups
    '''
    def __init__(self, C_CLASS, C_OPT):
        '''
        Initiation
        Inputs:
            C_OPT(CASE_OPT class)
        '''
        print(C_OPT)
        self.CASE = C_CLASS
        self.CASE_OPT = C_OPT
    
    def read_json_base(self, json_file):
        '''
        Read the base json file
        '''
        assert(os.access(json_file, os.R_OK))
        with open(json_file, 'r') as fin:
            self.base_options = json.load(fin)

    def create_group(self, features, base_dir, output_dir, base_name, bindings):
        '''
        create new group
        '''
        total = 1
        sub_totals = []
        sizes = []
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)  # make dir if not existing
        for feature in features:
            sub_totals.append(total)
            sizes.append(len(feature.get_values()))
            total *= len(feature.get_values())
        if bindings == []:
            # no bindings assigned, just follow the order in values
            for i in range(total):
                options = self.base_options.copy()
                options['name'] =  base_name
                for j in range(len(features)):
                    feature = features[j]
                    x_j = i//sub_totals[j] % len(feature.get_values())
                    values = feature.get_values()
                    options = Utilities.write_dict_recursive(options,\
                        feature.get_keys(), values[x_j])
                    if feature.if_abbrev_value():
                        name_appendix = get_name_appendix(feature.get_abbrev_value_options(), values[x_j]) # appendix by value
                    else:
                        name_appendix = feature.get_abbrev_strings()[x_j]  # appendix by string
                    options['name'] += ('_' + name_appendix)
                options["output directory"] = output_dir
                case_dir = create_case_with_json(options, self.CASE, self.CASE_OPT)
                json_path = os.path.join(case_dir, "case.json")
                with open(json_path, 'w') as fout:
                    json.dump(options, fout, indent=2)
        else:
            # follow the assigned binding
            for binding in bindings:
                options = self.base_options.copy()
                options['name'] =  base_name
                options["output directory"] = output_dir
                options['base directory'] = base_dir
                for j in range(len(features)):
                    feature = features[j]
                    x_j = binding[j]
                    values = feature.get_values()
                    options = Utilities.write_dict_recursive(options,\
                        feature.get_keys(), values[x_j])
                    if feature.if_abbrev_value():
                        name_appendix = get_name_appendix(feature.get_abbrev_value_options(), values[x_j]) # appendix by value
                    else:
                        name_appendix = feature.get_abbrev_strings()[x_j]  # appendix by string
                    options['name'] += ('_' + name_appendix)
                case_dir = create_case_with_json(options, self.CASE, self.CASE_OPT)
                json_path = os.path.join(case_dir, "case.json")
                with open(json_path, 'w') as fout:
                    json.dump(options, fout, indent=2)
            
        pass


def get_name_appendix(options, value):
    '''
    name appendix to a case
    options (list): key and scale
    value(int or float)
    '''
    key = options[0]
    scale = options[1]
    if type(value) == int:
        num_str = "%d" % (int(value * scale))
    else:
        num_str = "%.1f" % (value * scale)
    return key + num_str


def CreateGroup(json_path, CASE, CASE_OPT):
    '''
    create a new group
    Inputs:
        json_path -
    Returns:
        -
    '''
    group_opt = GROUP_OPT()
    group_opt.read_json(json_path)
    feature_opt = group_opt.values[1][0]
    group = GROUP(CASE, CASE_OPT)
    group.read_json_base(group_opt.get_base_json_path())
    if os.path.isdir(group_opt.get_output_dir()):
        is_pursue = input("Directory (%s) already exists, delete ? (y/n): " % group_opt.get_output_dir())
        if is_pursue == 'y':
            rmtree(group_opt.get_output_dir())
            os.mkdir(group_opt.get_output_dir())
        else:
            is_pursue = input("update ? (y/n)")
            if is_pursue != 'y':
                print("Aborting")
                return 0
    else:
        os.mkdir(group_opt.get_output_dir())
    copy(json_path, os.path.join(group_opt.get_output_dir(), "group.json"))
    group.create_group(*group_opt.to_create_group())


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
    elif _commend == 'create_group':
        CreateGroup(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()