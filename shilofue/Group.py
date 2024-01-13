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
    cases = FindCasesInDir(group_dir)
    Utilities.my_assert(len(cases)>0, FileExistsError, "%s doesn't have cases in it." % group_dir)
    case_list = [] 
    step_list = []
    time_list = []
    wallclock_list = []
    for case in cases:
        # pull out information
        log_file = os.path.join(case, 'output', 'log.txt')
        assert(os.path.isfile(log_file))
        try:
            last_step, last_time, last_wallclock = RunTimeInfo(log_file, quiet=True)
        except TypeError:
            last_step = -1
            last_time = -1.0
            last_wallclock = -1.0
        case_list.append(os.path.basename(case))
        step_list.append(int(last_step))
        time_list.append(float(last_time))
        wallclock_list.append(float(last_wallclock))
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
        self.steps = []
        self.times = []
        self.wallclocks = []
        self.ab_paths = []
        self.attrs = ['cases', 'steps', 'times', 'wallclocks']
        self.n_case = 0
        self.has_update = True
        # todo_diagram
        self.VISIT_OPTIONS = kwargs.get("VISIT_OPTIONS", None)

    def import_directory(self, _dir):
        '''
        Import from a directory, look for groups and cases
        Inputs:
            _dir (str): directory to import
        '''
        assert(os.path.isdir(_dir))
        case_list, step_list, time_list, wallclock_list = ReadBasicInfoGroup(_dir)
        self.n_case += len(case_list)
        self.cases += case_list
        self.steps += step_list
        self.times += time_list
        self.wallclocks += wallclock_list
        for _case in case_list:
            self.ab_paths.append(os.path.join(_dir, _case))
    
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

        odata = reader.export("", self.attrs, data_only=True)
        case_list = odata[:, 0]
        self.n_case = odata.shape[0]

        for i in range(len(case_list)):
            case = case_list[i]
            if case not in self.cases:
                self.import_data(odata[i, :])
        
        # set has_update to False
        self.has_update = False

    def import_data(self, o_data_1d):
        '''
        import data from 1d array
        o_data_1d (1d array)
        '''
        self.cases.append(o_data_1d[0])
        self.times.append(o_data_1d[1])
        self.steps.append(o_data_1d[2])
        self.wallclocks.append(o_data_1d[3])
        self.ab_paths.append('')

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
        for i in range(len(attrs_to_output)):
            _attr = attrs_to_output[i]
            header = headers[i]
            temp = np.array(getattr(self, _attr))
            if len(temp) > 0:
                # len(temp) == 0 marks an void field
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
        
    def write_file(self, o_path):
        '''
        write file
        Inputs:
            o_path (str): path of output
        '''
        attrs_to_output = ['cases', 'steps', 'times', 'wallclocks']
        headers = ['cases', 'steps', 'times (yr)', 'wallclocks (s)']

        with open(o_path, 'w') as fout: 
            self.write(fout, attrs_to_output, headers) 
        
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

    # error types
    class SummaryFileNotExistError(Exception):
        pass
         

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