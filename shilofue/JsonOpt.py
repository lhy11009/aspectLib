# -*- coding: utf-8 -*-
r"""This module defines, uses and documents options from a json file

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
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def Usage():
    print("\
This module defines, uses and documents options from a json file\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \
        ")


class JSON_OPT():
    """
    This class defines json options and offer interfaces to use
    and document them
    Args:
        keys (list of list of str): each member contains a list of keys from the top level
        descriptions (list of str): each member is a description of the option of the file
        types (list of type): each member is the type of the variable
        values (list of variable (type)): each member is the value of a variable
    """
    def __init__(self):
        """
        Initiation
        """
        pass

    def read_option(self, _path):
        """
        Read in json option from file path
        Inputs:
            _path (str): path of the json file
        """
        assert(os.access(_path, os.R_OK))
        # json.load()

    def add_key(self, description, _type, keys):
        """
        Add an option (key, describtion)
        Inputs:
            description (str)
            _type (type)
            keys (list): a list of keys from the top level
        """
        for key in keys:
            assert(type(key) == _type)
        self.keys.append(keys)
        self.descriptions.append(description)
        self.types.append(_type)
        pass

    def get_value(self, keys):
        """
        Get the value through a list of keys
        Inputs:
            keys(list): a list of keys from the top level
        """
        pass

    def document(self):
        """
        Print the documentation of this class.
        """
        # print()
        pass


def ShowAllOptions(_path):
    '''
    Show all options in a json file as list of keys
    This aims to aid the process of interpreting and coding
    the options of this kind of file.
    Inputs:
        _path (str): path of the input json file
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
    elif _commend == 'show_all_options':
        # example:
        ShowAllOptions(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()