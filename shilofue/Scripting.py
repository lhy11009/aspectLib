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
import sys, os, argparse, re
# import json
# import pathlib
# import subprocess
import numpy as np
import warnings
# from matplotlib import cm
from matplotlib import pyplot as plt
from shutil import copytree, rmtree

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
  - Generate script: \n\
\n\
        python -m shilofue.Scripting generate_script -i /home/lochy/ASPECT_PROJECT/aspectLib/shilofue/PlotStatistics.py -o .test/test_output_scripting\n\
\n\
        ")
    

class SCRIPTING():
    '''
    A class that deals with scripting
    '''
    def __init__(self, file_path):
        '''
        Initiation, read contents and header
        Attribute:
            c_header: header that deals with comments
            file_path: input file to convert
        '''
        self.c_header = []
        self.ex_header = []
        self.im_header = []
        self.contents = []
        fin = open(file_path, 'r')
        line = fin.readline()
        self.file_path = file_path
        self.c_header = ReadCommentingHeaders(file_path)
        # read the importing from an internal module
        self.im_header = ReadInHeaders1(file_path)
        # read the importing from an external module
        self.ex_header = ReadExHeaders(file_path)
        ex_additional = []  # additional commands
        ex_additional.append("current_dir = os.path.dirname(__file__)\n")
        ex_additional.append("JSON_FILE_DIR = os.path.join(current_dir, \'json_examples\')\n")
        ex_additional.append("sys.path.append(os.path.join(current_dir, \'utilities\', \'python_scripts\'))\n")
        ex_additional.append("import Utilities\n")
        self.ex_header += ex_additional
        # read the contents of the file
        self.contents = ReadContents(file_path)

    def __call__(self, o_path, **kwargs):
        '''
        Write to an output file
        '''
        parse_dependence = kwargs.get("parse_dependence", False)
        recursive = kwargs.get("recursive", False)

        # todo_script
        o_dir = os.path.dirname(o_path)
        o_path1 = os.path.join(o_dir, "dependence.py")
        fout = open(o_path, 'w')
        if parse_dependence:
            fout1 = open(o_path1, 'w')
        
        # write commenting header
        for header in self.c_header:
            fout.write(header)
        
        # write external header
        for header in self.ex_header:
            fout.write(header)
        
        # write internal header
        for header in self.im_header:
            # module, objects = ParseImportSyntax(header)
            with open(self.file_path, 'r') as fin:
                slines = fin.readlines()
            module, alias = ParseModuleObject(header) # parse the module and alias
            if recursive:
                # call the recursive import function
                modules = [module]
                aliases = [alias]
                all_objects = [[]]
                modules, aliases, all_objects = FindImportModuleRecursive(modules, aliases, all_objects, slines)
            else:
                # call the non-recursive import function
                objects = FindImportModule(module, alias, slines)
            # print("module:", module) # debug
            # print("alias:", alias)
            # print("objects:", objects)
            explicit_import_contents = ""
            for _object in objects:
                # write contents of explicit import
                explicit_import_contents += ExplicitImport(module, _object, alias=alias)
                if parse_dependence:
                    fout1.write(explicit_import_contents)
                else:
                    fout.write(explicit_import_contents)
            for _object in objects:
                # replace alias.object with alias_object
                for i in range(len(self.contents)):
                    self.contents[i] = re.sub("%s.%s" % (alias, _object), "%s_%s" % (alias, _object), self.contents[i])

        # write main script contents 
        for content in self.contents:
            fout.write(content)
        print("SCRIPTING: write scripting %s" % (o_path))  # screen output


####
# Utility functions
####
def ReadCommentingHeaders(file_path):
    '''
    Read the first part of the file header: comments of the scripts
    Inputs:
        file_path(str): path of the file
    Return:
        headers: headers that contain the comments
    '''
    headers = []
    fin = open(file_path, 'r')
    line = fin.readline()
    while line != "":
        if re.match("^from.*import", line) or re.match("^import", line)\
            or re.match("^def", line):
            # the next part of the file is hit: break
            break
        else:
            headers.append(line)
        line = fin.readline()
    return headers


def ReadInHeaders(file_path):
    '''
    Read the second part of the file header: importing module from inside the package
    Inputs:
        file_path(str): path of the file
    Return:
        headers: headers that contain the importing of internal modules
    '''
    headers = []
    fin = open(file_path, 'r')
    line = fin.readline()
    start = False
    while line != "":
        if re.match("^from.*import", line) or re.match("^import", line):
            start = True
        elif re.match("^def", line) or re.match("^class", line):
            # this indicates that the next part of the file is reached
            break
        if re.match("^from.*shilofue.*import", line):
            headers.append(re.sub('[ \t\n]*$', '', line))
        line = fin.readline()
    return headers


class HeaderIncludeError(Exception):
    pass


def ReadInHeaders1(file_path):
    '''
    Read the second part of the file header: importing module from inside the package
    Inputs:
        file_path(str): path of the file
    Return:
        headers: headers that contain the importing of internal modules
    '''
    headers = []
    fin = open(file_path, 'r')
    line = fin.readline()
    while line != "":
        if re.match("^import.*shilofue.*as", line):
            headers.append(re.sub('[ \t\n]*$', '', line))
        elif re.match("^from.*shilofue.*import.*", line) or re.match("^import.*shilofue", line):
            # containing not-allowed header
            warnings.warn("The shilofue modules need to be " + \
                          "imported by import shiofue.\{module\} as \{alias\}." + 
                          "The current line is :\n%s" % line)
        elif re.match("^def", line) or re.match("^class", line):
            # this indicates that the next part of the file is reached
            break
        line = fin.readline()
    return headers



def ReadExHeaders(file_path):
    '''
    Read the third part of the file header: importing module from outside the package
    Inputs:
        file_path(str): path of the file
    Return:
        headers: headers that contain the importing of external modules
    '''
    headers = []
    fin = open(file_path, 'r')
    line = fin.readline()
    start = False
    while line != "":
        if re.match("^from.*import", line) or re.match("^import", line):
            start = True
        elif re.match("^def", line) or re.match("^class", line):
            # this indicates that the next part of the file is reached
            break
        if (re.match("^from.*import", line) or re.match("^import", line)) and\
             (not re.match(".*shilofue.*", line)) and\
             (not re.match("^import.*Utilities", line)):
            # this indicates this is an importing command from an external module
            # the "from.*shilofue.*import" and "import.*Utilities" will be handled
            # separately
            headers.append(line)
        line = fin.readline()
    return headers


def ReadContents(file_path):
    '''
    Read the fourtch part of the file: the contents of the file
    Inputs:
        file_path(str): path of the file
    Return:
        contents: contents that contain the importing of external modules
    '''
    contents = []
    fin = open(file_path, 'r')
    line = fin.readline()
    start = False
    while line != "":
        if re.match("^def", line) or re.match("^class", line):
            start = True
        if start:
            contents.append(line)
        line = fin.readline()
    return contents


def ExplicitImport(module, _object, _type=None, **kwargs):
    '''
    explicitly import object from module
    Inputs:
        module: name of the module
        _object: name of the object
        _type: specify the type of the object (class / function)
    Returns:
        contents of the object
    '''
    alias = kwargs.get("alias", None)
    prefix = "(class|def)"
    if _type is not None:
        # specify the type of the object
        assert(_type in ["class", "function"])
        if _type == "class":
            prefix = "class"
        else:
            prefix = "def"
    # find the path of the file
    file_path = ASPECT_LAB_DIR
    for part in module.split('.'):
        file_path = os.path.join(file_path, part)
    file_path += '.py'
    Utilities.my_assert(os.path.isfile(file_path), FileExistsError, "%s doesn't exist" % file_path)
    # import the object
    fin = open(file_path, 'r')
    line = fin.readline()
    content_list = []
    while line != "":
        ctemp = "^%s %s" % (prefix, _object)
        # print("ctemp: ") # debug
        # print(ctemp)
        if re.match("^%s %s" % (prefix, _object), line):
            if alias is not None:
                # if alias is present, add it as a prefix to the function
                cline = line.replace(_object, "%s_%s" % (alias, _object))
                # cline = re.sub(_object, "%s_%s" % (alias, _object), line)
            else:
                cline = line
            content_list.append(cline)
            line = fin.readline()
            while line != "":
                if re.match("^(class|def)", line):
                    # next function is reached, stop
                    break
                content_list.append(line)
                line = fin.readline()
            break
        line = fin.readline()
    # remove the vacant lines at the end
    for i in range(len(content_list)-1, -1, -1):
        content = content_list[i]
        if re.match('(\t| )*\n$', content):
            content_list.pop(i)
        else:
            break
    # construct the contents
    contents = ""
    for content in content_list:
        contents += content
    return contents


def ParseImportSyntax(line):
    '''
    parse from an importing syntax
    Inputs:
        line: one line input
    Returns:
        module (str): name of the module
        _object (str): name of the object
    '''
    # line must have importing syntax
    assert(re.match("^from.*import", line))
    # parse module
    temp = re.sub("^from(\t| )*", '', line)
    temp = re.sub("(\t| )*import.*$", '', temp)
    module = re.sub("\n$", "", temp) # remove trailing newlines
    # parse object
    temp = re.sub("(\t| )*as.*$", '', line)  # remove the as syntax
    temp = re.sub(".*import(\t| )*", '', temp)
    objects = temp.split(',')
    for i in range(len(objects)):
        temp =  re.sub("(\t| )*", '', objects[i])
        objects[i] = re.sub("(\t| )*\n$", '', temp)
    print("module:", module)
    return module, objects


def FindImportModule(module, alias, slines):
    '''
    Find imported module based on what is included in the file
    Inputs:
        module (str): name of the module
        alias (str): alias of the module
        slines (list): contents to look for the module from
    Returns:
        objects (list): objects imported from the module in the
            current file contents
    '''
    assert(type(slines) == list)
    # then search for objects in the file
    objects = []
    for sline in slines:
        if re.match("^.*" + alias, sline) and (not re.match("import", sline)):
            temp = re.sub(".*" + alias + ".", "", sline)  # remove alias.
            temp = re.sub("\..*$", "", temp) # remove . (class functions)
            temp = re.sub("\).*$", "", temp) # remove right parenthesis (in case this is a class)
            temp = re.sub("\n", "", temp) # remove endline
            _object = re.sub("\(.*$", "", temp) # remove left parenthesis
            objects.append(_object)
    return objects


# todo_script
#def FindImportModuleRecursive(modules, aliases, all_objects, slines):
#    '''
#    Find imported module based on what is included in the file
#    Perform the lookup recursively until all the dependence is
#    found.
#    Inputs:
#        modules (str): name of the modules
#        aliases (str): aliases of the modules
#        all_objects (list): a list of all the objects imported in all modules
#        slines (list): contents to look for the module from
#    Returns:
#        modules (str): name of the modules
#        aliases (str): aliases of the modules
#        all_objects (list): a list of all the objects imported in all modules
#        slines (list): current file contents
#    '''
#    objects = FindImportModule(module, alias, slines)
#    explicit_import_contents = ""
#    for _object in objects:
#        # write contents of explicit import
#        explicit_import_contents += ExplicitImport(module, _object, alias=alias)
#    new_lines = explicit_import_contents.split('\n')
#    for header in headers
#        ParseModuleObject(headers)
#    for line in new_lines:
#
#    if  
#        return FindImportModuleRecursive(modules, aliases, all_objects, new_lines)
#    else:
#        return modules, aliases, all_objects


# todo_script
def ParseModuleObject(line):
    '''
    Parse module and object from a string input
    Inputs:
        line (str): string to parse the module and object
    Returns:
        module (str): name of the module
        alias (str): alias of the module
    '''
    # first read module from header
    assert(re.match("^import.*", line))
    temp = re.sub("^import(\t| )", '', line)
    module = re.sub("(\t| )+as.*$", '', temp)
    alias = re.sub(".*(\t| )+as(\t| )+", '', temp)
    return module, alias


def ParseHeader(file_path):
    '''
    Parse the header of the input file, figuring out
    what to import. Only deal with the "from ... import ..."
    syntax
    Inputs:
        file_path: path of the input file
    Returns:
        module_list: list of module
        object_list: list of object
    '''
    module_list = []
    object_list = []
    Utilities.my_assert(os.path.isfile(file_path), FileExistsError, "%s doesn't exist" % file_path)
    # read headers
    fin = open(file_path, 'r')
    line = fin.readline()
    headers = []
    while line != "":
        line = fin.readline()
        if re.match("^from", line):
            header = re.sub(' *(#.*)?\n$', '', line)
            headers.append(header)
        if re.match("^(class|def)", line):
            break
    # read module and object from header
    for header in headers:
        module, _object = ParseImportSyntax(header)
        module_list.append(module)
        object_list.append(_object)
    return module_list, object_list


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
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='Some outputs')
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
    elif _commend == 'generate_script':
        # inputs: a file
        # outputs: a directory that exists
        assert(os.path.isfile(arg.inputs))
        Scripting = SCRIPTING(arg.inputs)
        assert(os.path.isdir(arg.outputs))
        ofile = os.path.join(arg.outputs, os.path.basename(arg.inputs))
        Scripting(ofile, parse_dependence=True)
        assert(os.path.isfile(ofile))
        # copy the utilities
        utilities_dir = os.path.join(ASPECT_LAB_DIR, 'utilities')
        utilities_o_dir = os.path.join(arg.outputs, 'utilities')
        if os.path.isdir(utilities_o_dir):
            # remove old outputs
            rmtree(utilities_o_dir)
        assert(os.path.isdir(utilities_dir))
        copytree(utilities_dir, utilities_o_dir)
        assert(os.path.isdir(utilities_o_dir))
        # copy the json_files
        json_dir = os.path.join(ASPECT_LAB_DIR, 'files', "json_examples")
        json_o_dir = os.path.join(arg.outputs, 'json_examples')
        if os.path.isdir(json_o_dir):
            # remove old outputs
            rmtree(json_o_dir)
        assert(os.path.isdir(json_dir))
        copytree(json_dir, json_o_dir)
        assert(os.path.isdir(json_o_dir))

    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()