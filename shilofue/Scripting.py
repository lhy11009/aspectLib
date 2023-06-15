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

####
# Utility functions
####

def ExplicitImport(module, _object, _type=None):
    '''
    explicitly import object from module
    Inputs:
        module: name of the module
        _object: name of the object
        _type: specify the type of the object (class / function)
    Returns:
        contents of the object
    '''
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
        if re.match("^%s %s" % (prefix, _object), line):
            content_list.append(line)
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
    print(contents) # debug
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
    module = re.sub("(\t| )*import.*$", '', temp)
    # parse object
    temp = re.sub("(\t| )*as.*$", '', line)  # remove the as syntax
    _object = re.sub(".*import(\t| )*", '', temp)
    return module, _object


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
    elif _commend == 'foo':
        # example:
        SomeFunction('foo')
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()