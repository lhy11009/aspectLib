# -*- coding: utf-8 -*-
r"""Work with phase transitions in aspect

This exports: 

  -

This depends on:

  -  

""" 
import numpy as np
import sys, os, argparse
import json
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


class PHASE_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with phase transitions
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Reference density of this composition", float, ["rho0"], 3300.0, nick='rho0')
        self.add_key("Density differences", list, ["drho"], [], nick='drho')
        self.add_key("Compositional volume proportion", list, ["xc"], [], nick='xc')


class CDPT_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with composition-dependent phase transitions(CDPT)
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_features("Individual compositions", ["compositions"], PHASE_OPT, nick='compositions')


def ParsePhaseInput(inputs):
    '''
    parse input of phases to a aspect input form
    
    '''
    output = ""

    # read in density of the base phase
    rho_base = inputs.get("rho_base", 3300.0)
    Utilities.my_assert(type(rho_base) == float, TypeError, "base value for density must be a float")
    # read in the density change with each transtion
    drho = inputs.get("drho", [0.0])
    Utilities.my_assert(type(drho) == list, TypeError, "value of density change must be a list")
    # read in the fraction with each transtion
    xc = inputs.get("xc", [1.0])
    Utilities.my_assert(type(xc) == list, TypeError, "value of fraction with transition must be a list")
    # number of transition
    total = len(drho)
    Utilities.my_assert(len(xc) == total, ValueError, "length of xc and drho must be the same")

    # compute density
    rho = rho_base * np.ones(total + 1)
    for i in range(total):
        rho[i+1:] += drho[i] * xc[i]

    # generate output
    output += "%.1f" % rho[0]
    for i in range(1, total+1):
        output += "|%.1f" % rho[i]

    return output

def Usage():
    CDPT_opt=CDPT_OPT()
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
    - parse composition dependent phase transition inputs:\n\
\n\
        python -m shilofue.Parse phase_input -i ./files/TwoDSubduction/phases_1_0.json\n\
\n\
    %s\n\
\n\
" % CDPT_opt.document()
        )

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
    elif _commend == 'phase_input':
        # example:
        Utilities.my_assert(os.access(arg.inputs, os.R_OK), FileExistsError, "Json file doesn't exist.")
        with open(arg.inputs) as fin:
            inputs = json.load(fin)

        # get the outputs
        outputs = "density = "
        for key, value in inputs.items():
            if type(value) == dict:
                output = ParsePhaseInput(value)
                outputs += "%s: %s, " % (key, output)

        # print the output
        print(outputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()