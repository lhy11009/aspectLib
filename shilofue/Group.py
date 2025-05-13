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
from shilofue.PlotRunTime import RunTimeInfo
import shilofue.ParsePrm as ParsePrm
from shilofue.Plot import LINEARPLOT
# import pathlib
# import subprocess
import numpy as np
import pandas as pd
import subprocess
# from matplotlib import cm
from matplotlib import pyplot as plt
from shutil import rmtree, copy, SameFileError

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

HaMaGeoLib_DIR = "/home/lochy/ASPECT_PROJECT/HaMaGeoLib"
if os.path.abspath(HaMaGeoLib_DIR) not in sys.path:
    sys.path.append(os.path.abspath(HaMaGeoLib_DIR))
from hamageolib.utils.case_options import parse_log_file_for_time_info
from hamageolib.utils.plot_helper import fix_wallclock_time
from hamageolib.utils.file_reader import read_aspect_header_file

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \n\
\n\
  - Generate group documentation: \n\
    python -m shilofue.Group document -i /mnt/lochy0/ASPECT_DATA/ThDSubduction \n\
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
        self.add_key("if appending abbreviation (same length with values)", list, ["if abbreviating"], [], nick='if_append_abbrev')

    def check(self):
        '''
        check if everything is right
        '''
        _name = self.values[0]
        _values = self.values[2]
        abbrev_strings = self.values[5]
        assert(len(self.values[3]) == 2)  # abbreviation has 2 entries (name, scale)
        if self.if_abbrev_value() == 1:
            # contains a prefix and a scaling
            assert(len(self.values[3]) == 2)
        else:
            # number of strings equal number of values
            Utilities.my_assert(len(abbrev_strings) == len(_values), ValueError, "the options for abbrev_strings for \"%s\" is not correct." % _name)
        pass
        # check the values assigned for "if abbreviation"
        is_abbrev_value = self.values[4]
        assert(is_abbrev_value in [0, 1])
        if is_abbrev_value == 1:
            abbrev_value_options = self.values[3]  # abbrev_value_options contains a string and a float
            assert(len(abbrev_value_options)==2 and\
            type(abbrev_value_options[0]) == str and\
            type(abbrev_value_options[1]) == float)
        if self.values[6] == []:
            # no assigned, abbrevtion for all
            self.values[6] = [1 for i in range(len(self.values[2]))]
        else:
            Utilities.my_assert(len(self.values[6]) == len(self.values[2]), ValueError,\
            "length of the option (%s) is different from the other option (%s)" % (self.values[6], self.values[2]))
            for i in self.values[6]:
               Utilities.my_assert(i == 0 or i == 1, ValueError,\
               "Value of \"if abbreviation\" is either 0 or 1") 
            

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
    
    def if_append_abbrev(self, index):
        '''
        whether appending abbrevation to case name for entry i
        in the list of values
        Inputs:
            index (int): index of the entry
        '''
        return self.values[6][index]


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
        self.add_features('Base Feature to set up for the group', ['base features'], FEATURE_OPT)
        # todo_comb
        self.add_key('Combine case run', int, ['combine case run'], 0, nick='combine_case_run')
        self.add_key('Slurm options', list, ['slurm'], [], nick='slurm options')
    
    def check(self):
        if self.values[4] != []:
            for binding in self.values[4]:
                assert(len(binding) == len(self.values[1]))  # binding has length of features
                for i in range(len(binding)):
                    assert(type(binding[i]) == int)
        # check base features only have 1 value
        base_features = self.values[6]
        for base_feature in base_features:
            assert(len(base_feature.get_values())==1)

        # check the slurm file paths
        slurm_options = self.values[8]
        for slurm_option in slurm_options:
            slurm_file_base =  Utilities.var_subs(slurm_option["slurm file"])
            assert(os.path.isfile(slurm_file_base))

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
        base_features = self.values[6]
        features = self.values[1]
        combine_case_run = self.values[7]
        slurm_options = self.values[8]
        return base_features, features, Utilities.var_subs(self.values[5]),\
        Utilities.var_subs(self.values[3]), self.values[0], self.values[4],\
        combine_case_run, slurm_options


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
        Utilities.my_assert(os.access(json_file, os.R_OK), FileExistsError, "%s doesn't exist." % json_file)
        with open(json_file, 'r') as fin:
            self.base_options = json.load(fin)

    def create_group(self, base_features, features, base_dir, output_dir, base_name, bindings,\
                     combine_case_run, slurm_options, **kwargs):
        '''
        create new group
        Inputs:
            kwargs (dict):
                update (bool): update existing cases?
        '''
        total = 1
        sub_totals = []
        sizes = []
        is_update = kwargs.get('update', False)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)  # make dir if not existing
        # handle base features, set up group-wise value of varibles
        options = self.base_options.copy()
        for base_feature in base_features:
            options = Utilities.write_dict_recursive(options,\
                    base_feature.get_keys(), base_feature.get_values()[0])
        self.base_options = options
        # todo_comb
        # slurm configurations
        self.base_options['slurm'] = slurm_options
        
        # handle features, varying the value of variables in cases
        case_names = []
        for feature in features:
            sub_totals.append(total)
            sizes.append(len(feature.get_values()))
            total *= len(feature.get_values())
        if bindings == []:
            # no bindings assigned, just follow the order in values
            for i in range(total):
                options = self.base_options.copy()
                options['name'] =  base_name
                # assemble the features for case option
                for j in range(len(features)):
                    feature = features[j]
                    x_j = i//sub_totals[j] % len(feature.get_values())
                    values = feature.get_values()
                    options = Utilities.write_dict_recursive(options,\
                        feature.get_keys(), values[x_j])
                    if feature.if_append_abbrev(x_j):
                        if feature.if_abbrev_value():
                            name_appendix = get_name_appendix(feature.get_abbrev_value_options(), values[x_j]) # appendix by value
                        else:
                            name_appendix = feature.get_abbrev_strings()[x_j]  # appendix by string
                        options['name'] += ('_' + name_appendix)
                options["output directory"] = output_dir
                # call function to create the case
                case_dir = create_case_with_json(options, self.CASE, self.CASE_OPT, update=is_update, fix_base_dir=base_dir)
                case_names.append(os.path.basename(case_dir))
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
                # assemble the features for case option
                for j in range(len(features)):
                    feature = features[j]
                    x_j = binding[j]
                    values = feature.get_values()
                    try:
                        options = Utilities.write_dict_recursive(options,\
                            feature.get_keys(), values[x_j])
                    except IndexError:
                        raise IndexError("values have length %d, while index is %d\
It's likely a wrong value is assigned in the \"binding\" part of the json file." % (len(values), x_j))
                    if feature.if_append_abbrev(x_j):
                        if feature.if_abbrev_value():
                            name_appendix = get_name_appendix(feature.get_abbrev_value_options(), values[x_j]) # appendix by value
                        else:
                            name_appendix = feature.get_abbrev_strings()[x_j]  # appendix by string
                        options['name'] += ('_' + name_appendix)
                # call function to create the case
                case_dir = create_case_with_json(options, self.CASE, self.CASE_OPT, update=is_update)
                case_names.append(os.path.basename(case_dir))
                json_path = os.path.join(case_dir, "case.json")
                with open(json_path, 'w') as fout:
                    json.dump(options, fout, indent=2)

        # todo_comb 
        # Combine case run in slurm
        if combine_case_run:
            for slurm_option in slurm_options:
                GenerateSlurmCombine(slurm_option, output_dir, case_names)



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
    elif value > 1.0:
        num_str = "%.1f" % (value * scale)
    elif value < 1.0:
        num_str = "%.2e" % (value * scale)
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
    is_update_existing = False
    if os.path.isdir(group_opt.get_output_dir()):
        is_pursue = input("Directory (%s) already exists, delete ? (y/n): " % group_opt.get_output_dir())
        if is_pursue == 'y':
            is_pursue1 = input("Delete, are you sure ? (y/n): ")
            if is_pursue1 == 'y':
                rmtree(group_opt.get_output_dir())
                os.mkdir(group_opt.get_output_dir())
        else:
            is_pursue = input("update ? (y/n)")
            if is_pursue != 'y':
                print("Aborting")
                return 0
            temp = input("update existing cases ? (y/n)")
            if temp == 'y':
                is_update_existing = True
            else:
                is_update_existing = False
    else:
        os.mkdir(group_opt.get_output_dir())
    try:
        copy(json_path, os.path.join(group_opt.get_output_dir(), "group.json"))
    except SameFileError:
        print("Copy: these are the same file (%s, %s), pass."\
        % (json_path, os.path.join(group_opt.get_output_dir(), "group.json")))
    
    group.create_group(*group_opt.to_create_group(), update=is_update_existing)


####
# classes and functions for slurm configuration
# todo_comb
####
def GenerateSlurmCombine(slurm_option, _dir, case_names, **kwargs):
    '''
    Inputs:
        slurm_options - options for slurm file
        _dir - output directory (the group directory)
        cases_names - names of cases in a group
        wargs:
            slurm_case_name: name of the slurm case
    ''' 
    size = len(case_names)
    slurm_file_base = slurm_option["slurm file"]
    nproc = slurm_option["cpus"]
    build_dir = slurm_option["build directory"]
    aspect_exec = "${ASPECT_SOURCE_DIR}/build_%s/aspect" % build_dir
    slurm_case_name = kwargs.get("slurm_case_name", os.path.basename(_dir))  # take the group name
    slurm_file_base_base = os.path.basename(slurm_file_base)

    # generate the slurm file
    slurm_file_path = os.path.join(_dir, slurm_file_base_base)
    SlurmOperator = ParsePrm.SLURM_OPERATOR(slurm_file_base)
    SlurmOperator.SetAffinity(1, nproc * size, 1)
    SlurmOperator.ResetCommand(command_only=1)
    SlurmOperator.SetName(slurm_case_name)
    SlurmOperator.SetTimeByHour(300)

    # generate the command to run
    extra_contents = ""
    # directories to run
    temp = "subdirs=("
    i = 0
    for case_name in case_names:
        if i > 0:
            temp += " "
        temp += "\"%s\"" % case_name
        i += 1
    temp += ")\n"
    extra_contents += temp
    # command
    exec =  "srun --exclusive --ntasks %d %s case.prm &" % (nproc, aspect_exec)
    # generate the loop
    temp = """
for subdir in ${subdirs[@]}; do
    cd ${subdir}
    %s
    cd ..
done
wait
""" % (exec)
    extra_contents += temp
    SlurmOperator.SetExtra(extra_contents)

    # generate file
    SlurmOperator(slurm_file_path)


####
# classes and functions for documenting
####
class GDOC_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with DOC
    List of keys:
    '''
    def __init__(self):
        self.add_key("Base directory to document", str, ["base directory"], ".", nick='base_dir')
        pass

    def to_execute(self):
        base_dir = self.values["base_dir"]
        return base_dir



class GDOC():
    '''
    generate documentation
    '''
    def __init__(self):
        '''
        Initiation
        '''
        self.groups = []
        self.group_names = []
        self.cases = []
        self.case_names = []
        self.steps = []
        self.times = []
        self.wallclocks = []

    def execute(self, base_dir, **kwargs):
        '''
        execute the documentation
        '''
        o_dir = kwargs.get('o_dir', base_dir)
        # create a documentation directory
        doc_dir = os.path.join(o_dir, "documentation")
        if not os.path.isdir(doc_dir):
            os.mkdir(doc_dir)
        # find all groups and seperate cases
        self.groups = Utilities.SortByCreationTime(FindGroupsInDir(base_dir))
        self.group_names = [os.path.basename(_path) for _path in self.groups]  # group names
        self.cases = Utilities.SortByCreationTime(FindCasesInDir(base_dir))
        self.case_names = [os.path.basename(_path) for _path in self.cases]  # group names
        # generate markdown files
        outputs = self.create_markdown()
        o_path = os.path.join(doc_dir, "group_doc.mkd")
        with open(o_path, 'w') as fout:
            fout.write(outputs)
        print("File %s generated." % o_path)
        o_path_html = os.path.join(doc_dir, "group_doc.html")
        # convert to html
        subprocess.run("pandoc %s -t html -o %s" % (o_path, o_path_html), shell=True)
        assert(os.path.isfile(o_path_html))
        print("Html file %s generated" % o_path_html)
        # generate latex files
        outputs = self.create_latex()
        o_path = os.path.join(doc_dir, "group_doc.tex")
        with open(o_path, 'w') as fout:
            fout.write(outputs)
        print("File %s generated." % o_path)

    def create_markdown(self):
        """
        Create contents of a markdown file
        """
        outputs = ""
        for i in range(len(self.groups)):
            group = self.groups[i]
            group_name = self.group_names[i]
            outputs += "### %s\n\n" % group_name
            outputs += self.output_case_table_tex(group) + "\n"
        return outputs
        
    def create_latex(self):
        '''
        Create contents of a latex file
        '''
        outputs = ""
        for i in range(len(self.groups)):
            group = self.groups[i]
            group_name = self.group_names[i]
            # outputs += "### %s\n\n" % group_name
            outputs += self.output_case_table_tex(group, format='latex') + "\n"
        return outputs


    class CreateCaseTableError(Exception):
        pass


    def output_case_table(self):
        '''
        case table output
        '''
        pass


    def output_case_table_tex(self, group_dir, **kwargs):
        '''
        create a table for case documentation
        Inputs:
            group_dir: directory of the group
            kwargs:
                format - output format
        Returns:
            table_contnents: contents of the table
        '''
        json_path = os.path.join(group_dir, 'group.json')
        _format = kwargs.get("format", "markdown")
        # read group options from the json file
        group_opt = GROUP_OPT()
        group_opt.read_json(json_path)
        features = group_opt.get_features()
        n_features = len(features)
        # find all cases
        cases = FindCasesInDir(group_dir)
        Utilities.my_assert(len(cases)>0, FileExistsError, "%s doesn't have cases in it." % group_dir)
        header = ["cases", "steps", "time", "wallclock"]
        colors = [None, "blue", "blue", "blue"]
        # colors =  
        # construct data
        data = []
        case_list, step_list, time_list, wallclock_list = ReadBasicInfoGroup(group_dir)
        data.append(case_list)
        data.append(step_list)
        data.append(time_list)
        data.append(wallclock_list)
        TexTable = Utilities.TEX_TABLE("table-%s" % os.path.basename(group_dir),\
                                        header=header, data=data, colors=colors) # class initiation
        table_contents = TexTable(format=_format)
        return table_contents


def ReadBasicInfoGroup(group_dir):
    '''
    Read basic information from a group
    Inputs:
        group_dir(str): directory contains a few cases
    Return:
        case_list, step_list, time_list, wallclock_list (list): list of outputs
    '''
    print("ReadBasicInfoGroup: reading group %s" % group_dir)
    cases = FindCasesInDir(group_dir)
    Utilities.my_assert(len(cases)>0, FileExistsError, "%s doesn't have cases in it." % group_dir)
    case_list = [] 
    step_list = []
    time_list = []
    wallclock_list = []
    for case in cases:
        # pull out information
        last_step = -1
        last_time = -1.0
        last_wallclock = -1.0
        log_file = os.path.join(case, 'output', 'log.txt')
        if os.path.isfile(log_file):
            try:
                last_step, last_time, last_wallclock = RunTimeInfo(log_file, quiet=True)
            except TypeError:
                pass
        case_list.append(os.path.basename(case))
        step_list.append(int(last_step))
        time_list.append(float(last_time))
        wallclock_list.append(float(last_wallclock))
        print(" Found case %s" % os.path.basename(case))
    return case_list, step_list, time_list, wallclock_list


def DocumentGroupsInDir(_dir):
    '''
    generate documentation for groups in a directory
    '''
    assert(os.path.isdir(_dir))
    GDoc = GDOC()
    GDoc.execute(_dir)
    pass


class CASE_SUMMARY():
    '''
    Attributes:
        n_case (int): number of cases
        cases: name of cases
        steps: end steps of cases
        times: end times of cases
        wallclocks: running time of cases on the wall clock
        ab_paths: absolution_paths of cases
        has_update: has update to perform
        VISIT_OPTIONS: the VISIT_OPTIONS class. This could be initiated later in the
            scope to parse case informations.
    '''
    def __init__(self, **kwargs):
        '''
        Initiation
        Inputs:
            kwargs:
                VISIT_OPTIONS - the VISIT_OPTIONS class
        '''
        self.cases = []
        self.names = []
        self.steps = []
        self.times = []
        self.wallclocks = []
        self.ab_paths = []
        self.includes = []
        self.attrs = ['includes', 'cases', 'names', 'steps', 'times', 'wallclocks']
        self.n_case = 0

        self.attrs_to_output = ['names', 'includes', 'steps', 'times', 'wallclocks']
        self.headers = ['names', 'includes', 'steps', 'times (yr)', 'wallclocks (s)']

        self.has_update = True
        self.VISIT_OPTIONS = kwargs.get("VISIT_OPTIONS", None)

    def export(self, _name):
        return np.array([float(i) for i in getattr(self, _name)])

    def __call__(self, _dir, **kwargs):
        '''
        Inputs:
            _dir (str): directory to import
        '''
        format = kwargs.get("format", "txt")
        # try to import a previous summary first
        # if not found, import the cases in the directory
        if format == "txt":
            i_name = 'case_summary.txt'
        elif format == "csv":
            i_name = 'case_summary.csv'
        else:
            raise NotImplementedError()
        i_path = os.path.join(_dir, i_name)

        if os.path.isfile(i_path):
            if format == "txt":
                self.import_txt(i_path)
            elif format == "csv":
                self.import_csv(i_path)
            self.Update(**kwargs)
        else:
            self.import_directory(_dir, **kwargs)

    def Update(self, **kwargs):
        '''
        Update on properties
        Inputs:
            kwargs
        '''
        hr = 3600.0 # s in hr
        # append the include field
        if len(self.includes) < self.n_case:
            self.includes += [1 for i in range(self.n_case - len(self.includes))]
        # update step, time and wallclock time
        for i in range(self.n_case):
            # todo_time
            step, _time, wallclock = -1, -1.0, -1.0
            # log_file = os.path.join(self.ab_paths[i], 'output', 'log.txt')
            # if os.path.isfile(log_file):
            #     try:
            #         step, _time, wallclock = RunTimeInfo(log_file, quiet=True)
            #     except TypeError:
            #         pass
            
            log_path = os.path.join(self.ab_paths[i], 'output', 'log.txt')
            
            if os.path.isfile(log_path):

                output_dir = os.path.join(self.ab_paths[i], "awk_outputs")
                if not os.path.isdir(output_dir):
                    os.mkdir(output_dir)

                output_path = os.path.join(output_dir, "run_time_output.txt") 

                try:
                    # parse time info
                    parse_log_file_for_time_info(log_path, output_path)

                    pd_data = read_aspect_header_file(output_path)
                except TypeError:
                    pass
                else:
                    # fix wallclock time
                    pd_data = fix_wallclock_time(pd_data)
                    step = pd_data["Time step number"].iloc[-1]
                    _time = pd_data["Time"].iloc[-1]
                    wallclock = pd_data["Wall Clock"].iloc[-1]

            self.steps[i] = step
            if _time > 0.0:
                self.times[i] = _time / 1e6
            else:
                self.times[i] = -1.0
            if wallclock > 0.0:
                self.wallclocks[i] = wallclock / hr
            else:
                self.wallclocks[i] = -1.0
        # These fields need to be mannualy assigned, so we
        # only initiation a nan value for the first time
        if self.names == []:
            self.names = [np.nan for i in range(self.n_case)]
        else:
            pass

    def import_directory(self, _dir, **kwargs):
        '''
        Import from a directory, look for groups and cases
        Inputs:
            _dir (str): directory to import
        '''
        hr = 3600.0 # s in hr
        assert(os.path.isdir(_dir))
        # first parse in run time information
        case_list, step_list, time_list, wallclock_list = ReadBasicInfoGroup(_dir)
        for i in range(len(case_list)):
            _case = case_list[i]
            step = step_list[i]
            _time = time_list[i]
            wallclock = wallclock_list[i]
            if _case in self.cases:
                j = self.cases.index(_case)
                self.steps[j] = step
                if _time > 0.0:
                    self.times[j] = _time / 1e6
                else:
                    self.times[j] = -1.0
                if wallclock > 0.0:
                    self.wallclocks[j] = wallclock / hr
                else:
                    self.wallclocks[j] = -1.0
            else:
                self.n_case += 1
                self.cases.append(_case)
                self.steps.append(step)
                if _time > 0.0:
                    self.times.append(_time / 1e6)
                else:
                    self.times.append(-1.0)
                if wallclock > 0.0:
                    self.wallclocks.append(wallclock / hr)
                else:
                    self.wallclocks.append(-1.0)
                self.ab_paths.append(os.path.join(_dir, _case))
        
        # append the include field
        if len(self.includes) < self.n_case:
            self.includes += [1 for i in range(self.n_case - len(self.includes))]

    def import_txt(self, i_path):
        '''
        import an existing txt file
        Inputs:
            i_path: path to a txt file
        '''
        reader = LINEARPLOT('case_summary')

        Utilities.my_assert(os.path.isfile(i_path), self.SummaryFileNotExistError, "%s is not a file." % i_path)
        reader.ReadHeader(i_path)
        reader.ReadData(i_path, dtype=str)

        # Read input cases
        # if case already exists in the record, update
        # else: append new case into the record
        n_case_old = self.n_case
        in_cases = reader.export_field_as_array("cases")
        in_ab_paths = reader.export_field_as_array("ab_paths")
        indexes = []
        for i in range(len(in_cases)):
            in_case = in_cases[i]
            if in_case in self.cases:
                j = self.cases.index(in_case)
                indexes.append(j)
            else:
                self.cases.append(in_case)
                self.ab_paths.append(in_ab_paths[i])
                indexes.append(-1)
                self.n_case += 1

        # import data and either append or update fields
        for attr in self.attrs:
            # "cases" attribute is handled in the previous loop
            if attr == "cases":
                continue
            temp = getattr(self, attr)
            in_array = reader.export_field_as_array(attr)
            if (len(in_array) == 0):
                # no data input, continue
                continue
            if len(temp) > 0:
                # field exist
                for i in range(len(in_cases)):
                    if indexes[i] < 0:
                        # case doesn't preexist
                        temp.append(in_array[i])
                    else:
                        # case preexist
                        temp[indexes[i]] = in_array[i] 
            else:
                # field doesn't exist
                # current: do nothing
                # TODO: figure out updating strategy
                for i in range(len(in_cases)):
                    if indexes[i] < 0:
                        # case doesn't preexist
                        # compensate results for previous case by adding trivial values
                        # then append the values for the new cases
                        if i == 0:
                            temp = [-1 for j in range(n_case_old)]
                        temp.append(in_array[i])
                    else:
                        # case preexist
                        # temp[indexes[i]] =in_array[i] 
                        continue
                # field doesn't exist
            setattr(self, attr, temp)
        
        # set has_update to False
        self.has_update = False

    def import_file(self, i_path):
        '''
        import an existing file
        Inputs:
            i_path: path to a txt/csv file
        '''
        reader = LINEARPLOT('case_summary') # initiation, reading txt
        df = None # initiation, pd object

        Utilities.my_assert(os.path.isfile(i_path), self.SummaryFileNotExistError, "%s is not a file." % i_path)

        format = i_path.split('.')[-1]
        if format == 'txt':
            reader.ReadHeader(i_path)
            reader.ReadData(i_path, dtype=str)
        elif format == 'csv':
            df = pd.read_csv(i_path, dtype=str)
        else:
            raise TypeError("%s: input file needs to be \"txt\" or \"csv\", but filename is %s" % Utilities.func_name(), i_path)

        # Read input cases
        # if case already exists in the record, update
        # else: append new case into the record
        n_case_old = self.n_case
        in_cases = None; in_ab_paths = None
        if format == 'txt':
            in_cases = reader.export_field_as_array("cases")
            in_ab_paths = reader.export_field_as_array("ab_paths")
        elif format == 'csv':
            in_cases = df["cases"]
            in_ab_paths = df["ab_paths"]
        indexes = []
        for i in range(len(in_cases)):
            in_case = in_cases[i]
            if in_case in self.cases:
                j = self.cases.index(in_case)
                indexes.append(j)
            else:
                self.cases.append(in_case)
                self.ab_paths.append(in_ab_paths[i])
                indexes.append(-1)
                self.n_case += 1

        # import data and either append or update fields
        for attr in self.attrs:
            # "cases" attribute is handled in the previous loop
            if attr == "cases":
                continue
            temp = getattr(self, attr)
            in_array = None
            if format == 'txt':
                in_array = reader.export_field_as_array(attr)
            if format == 'csv':
                in_array = None
                try:
                    in_array = df[attr]
                except KeyError:
                    in_array = []
            if (len(in_array) == 0):
                # no data input, continue
                continue
            if len(temp) > 0:
                # field exist
                for i in range(len(in_cases)):
                    if indexes[i] < 0:
                        # case doesn't preexist
                        temp.append(in_array[i])
                    else:
                        # case preexist
                        temp[indexes[i]] = in_array[i] 
            else:
                # field doesn't exist
                # current: do nothing
                # TODO: figure out updating strategy
                for i in range(len(in_cases)):
                    if indexes[i] < 0:
                        # case doesn't preexist
                        # compensate results for previous case by adding trivial values
                        # then append the values for the new cases
                        if i == 0:
                            temp = [-1 for j in range(n_case_old)]
                        temp.append(in_array[i])
                    else:
                        # case preexist
                        # temp[indexes[i]] =in_array[i] 
                        continue
                # field doesn't exist
            setattr(self, attr, temp)
        
        # set has_update to False
        self.has_update = False

        print("%s: file %s imported successfully" % (Utilities.func_name(), i_path))
    
    def write(self, fout, attrs_to_output, headers):
        '''
        Write an output file
        Inputs:
            fout (object): object of output
            attrs_to_output: attribute to output
            headers: headers of the file
        '''
        # header and data to write
        data_raw = []
        oheaders = []
        length_of_data = 0
        for i in range(len(attrs_to_output)):
            _attr = attrs_to_output[i]
            header = headers[i]
            temp = np.array(getattr(self, _attr))
            if i == 0:
                length_of_data = len(temp)
            if len(temp) > 0:
                # len(temp) == 0 marks an void field
                Utilities.my_assert(length_of_data == len(temp), ValueError,\
                    "Data for field \'%s\'(%d) doesn't much the length of data (%d)" % (_attr, len(temp), length_of_data))
                data_raw.append(temp)
                oheaders.append(header)
        data = np.column_stack(data_raw)
        # output
        # write header
        i = 0
        for header in oheaders:
            fout.write("# %d: %s\n" % (i+1, header))
            i += 1
        np.savetxt(fout, data, delimiter=" ", fmt="%s")

    def write_csv(self, fout, attrs_to_output, headers, mask=None):
        '''
        Write an output file in csv format
        Inputs:
            fout (object): object of output
            attrs_to_output: attribute to output
            headers: headers of the file
        '''
        # header and data to write
        data = {} # a vacant dict
        length_of_data = 0
        for i in range(len(attrs_to_output)):
            _attr = attrs_to_output[i]
            header = headers[i]
            temp = np.array(getattr(self, _attr))
            if i == 0:
                length_of_data = len(temp)
            if len(temp) > 0:
                # len(temp) == 0 marks an void field
                Utilities.my_assert(length_of_data == len(temp), ValueError,\
                    "Data for field \'%s\'(%d) doesn't much the length of data (%d), " % (_attr, len(temp), length_of_data)\
                        + str(temp))
                if mask is not None:
                    data[header] = temp[mask]
                else:
                    data[header] = temp
        df = pd.DataFrame(data)
        df.to_csv(fout) # save to csv file

    def sort_out_case_attribute_by_absolution_path(self, full_path, _attr):
        '''
        Sort out the case name assigned to a case path
        '''
        value = None
        for i in range(self.n_case):
            if self.ab_paths[i] == full_path:
                value = getattr(self, _attr)[i]
        return value

    def sort_out_by_attr_value(self, o_path, _attr, value):
        '''
        Sort out cases that have some value for an attribute
        Inputs:
            o_path (str) - file to output
            _attr (str or list) - name of the attribute
            value (float or list) - value of the attribute
        '''
        # append the case name first and the absolute path last
        attrs_to_output = ["cases"] + self.attrs_to_output + ["ab_paths"]
        headers = ["cases"] + self.attrs_to_output + ["ab_paths"]

        # handle two different situation:
        # 1. give 1 value
        # 2. give a list of values
        mask = None
        if type(_attr) == str: 
            temp = np.array(getattr(self, _attr))
            mask = np.abs((temp - value) / value) < 1e-6
        if type(_attr) == list:
            assert(len(_attr) == len(value))
            for i in range(len(_attr)):
                a = _attr[i]
                v = value[i]
                temp = np.array(getattr(self, a))
                mask1 = np.abs((temp - v) / v) < 1e-6
                if i == 0:
                    mask = mask1
                else:
                    mask = (mask & mask1)
        
        # write file
        format = o_path.split('.')[-1]
        if format == "csv":
            with open(o_path, 'w') as fout: 
                self.write_csv(fout, attrs_to_output, headers, mask)
        else:
            raise NotImplementedError()
        Utilities.my_assert(os.path.isfile(o_path), FileNotFoundError, "%s: file %s is not written successfully" % (Utilities.func_name(), o_path)) 
        print("%s: Write file %s" % (Utilities.func_name(), o_path))
        
    def write_file(self, o_path):
        '''
        write file
        Inputs:
            o_path (str): path of output
        '''
        # append the case name first and the absolute path last
        attrs_to_output = ["cases"] + self.attrs_to_output + ["ab_paths"]
        headers = ["cases"] + self.attrs_to_output + ["ab_paths"]

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

    def write_file_if_update(self, o_path):
        '''
        write file if there are updates
        Inputs:
            o_path (str): path of output
        '''
        if self.has_update:
            self.write_file(o_path)
        else:
            print("%s: File already exists %s" % (Utilities.func_name(), o_path))

    # todo_export 
    def export_to_latex_table(self, fout, attrs_to_output, headers):
        '''
        export summary results to a table in latex format
        '''
        # header and data to write
        data_raw = []
        oheaders = []
        length_of_data = 0
        for i in range(len(attrs_to_output)):
            _attr = attrs_to_output[i]
            header = headers[i]
            temp = np.array(getattr(self, _attr))
            if i == 0:
                length_of_data = len(temp)
            if len(temp) > 0:
                # len(temp) == 0 marks an void field
                Utilities.my_assert(length_of_data == len(temp), ValueError,\
                    "Data for field \'%s\'(%d) doesn't much the length of data (%d)" % (_attr, len(temp), length_of_data))
                data_raw.append(temp)
                oheaders.append(header)
        data = np.column_stack(data_raw)
        TexTable = Utilities.TEX_TABLE("case-summary", header=oheaders, data=np.transpose(data).tolist()) # class initiation
        table_contents = TexTable(format="latex", fix_underscore_in_content=False)
        fout.write(table_contents)

    # error types
    class SummaryFileNotExistError(Exception):
        pass

    def mesh_data(self, x_name, y_name, z_name, **kwargs):
        '''
        plot diagram
        Inputs:
            x_name, y_name, z_name - name of fields
            kwargs:
                nx: number of point along x
                ny: number of point along y
        Returns:
            XX, YY, ZZ - meshed data
        '''
        nx = kwargs.get("nx", 20)
        ny = kwargs.get("ny", 30)
        
        # get the properties
        Xs = getattr(self, x_name)
        Ys = getattr(self, y_name)
        Zs = getattr(self, z_name)
        x_min = np.min(np.array(Xs))
        x_max = np.max(np.array(Xs))
        y_min = np.min(np.array(Ys))
        y_max = np.max(np.array(Ys))

        # map to a 2d mesh
        x_temp = np.linspace(x_min, x_max, nx) 
        y_temp = np.linspace(y_min, y_max, ny)
        XX, YY = np.meshgrid(x_temp, y_temp)
        ZZ = WeightedByDistance(Xs, Ys, Zs, XX, YY)

        return XX, YY, ZZ
    
    def update_geometry(self, i):
        '''
        update info of geometry
        '''
        case_dir = self.ab_paths[i]
        Utilities.my_assert(os.path.isdir(case_dir), FileExistsError, "Directory doesn't exist %s" % case_dir)

        try:
            Visit_Options = self.VISIT_OPTIONS(case_dir)
            Visit_Options.Interpret()
            self.geometries[i] = Visit_Options.options["GEOMETRY"]
        except FileNotFoundError:
            self.geometries[i] = -1.0


def WeightedByDistance(Xs, Ys, Zs, XX, YY):
    '''
    Inputs:
        Xs, Ys - coordinates arrays
        Zs - data array
        XX, YY - meshed coordinates
    Return:
        zz - meshed data
    '''
    # initiate
    data_size = len(Xs)
    assert(data_size == len(Ys) and data_size == len(Zs))
    assert(XX.shape == YY.shape)
    ZZ = np.zeros(XX.shape)
    weight = np.zeros(XX.shape)

    # distance weight the data
    for i in range(data_size):
        dist = ((XX - Xs[i])**2.0 + (YY - Ys[i])**2.0)**0.5
        weight += 1.0 / dist
        ZZ += Zs[i] / dist
    ZZ /= weight
    return ZZ
         

####
# Utility fucntions
####
def FindCasesInDir(dir_path):
    '''
    Find cases within a directory
    '''
    cases = []
    for _name in os.listdir(dir_path):
        sub_dir_path = os.path.join(dir_path, _name)
        if os.path.isdir(sub_dir_path):
            prm_path = os.path.join(sub_dir_path, "case.prm")
            json_path = os.path.join(sub_dir_path, "case.json")
            if os.path.isfile(prm_path) and os.path.isfile(json_path):
                cases.append(sub_dir_path)
    return cases


def FindGroupsInDir(dir_path):
    '''
    Find groups within a directory
    '''
    groups = []
    for subdir, dirs, _ in os.walk(dir_path):
        for _dir in dirs:
            dir_path = os.path.join(subdir, _dir)
            json_path = os.path.join(dir_path, "group.json")
            if os.path.isfile(json_path):
                groups.append(dir_path)
    return groups



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
    elif _commend == 'document':
        DocumentGroupsInDir(arg.inputs)
        pass
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()