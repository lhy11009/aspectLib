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
        self.add_key("Name of this composition", str, ["name"], 'pyrolite', nick='name')
    
    def get_name(self):
        '''
        get name
        return: 
            name (str)
        '''
        return self.values[3]
    
    def get_dict(self):
        '''
        Return variables of PZ in a dict
        return:
            _dict (dict)
        '''
        _dict = {
            "rho0": self.values[0],
            "drho": self.values[1],
            "xc": self.values[2],
            "name": self.values[3]
        }
        return _dict
    


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
        # todo
        self.add_key("Apply and effective phase transition on the 660 for latent heat computation",\
        int, ["Effective 660 for latent heat"], 0, nick="effect_660_lh")
    
    def check(self):
        '''
        Check variables
        '''
        Utilities.my_assert(len(self.values[0]) > 0, ValueError, "Read parameters from file failed")
    
    def toPPTF(self):
        '''
        An interface to the ParsePhaseTransitionFile
        '''
        return self.values[0]


def ParsePhaseInput(inputs):
    '''
    parse input of phases to a aspect input form
    
    '''
    output = ""

    # number of transition
    total = len(inputs['drho'])
    rho = Get_reference_density(inputs)
    # generate output
    output += "%.1f" % rho[0]
    for i in range(1, total+1):
        output += "|%.1f" % rho[i]

    return output
        
def ParsePhaseTransitionFile(inputs):
    '''
    Parse the inputs of phase transition to output in aspect
    Inputs:
        inputs (dict): dictionary of phase transitions
    Returns:
        outputs (str): string for a aspect prm file
    '''
    outputs = "" 
    if type(inputs) == list:
        # read dict
        # get the outputs
        outputs = "density = "
        for phase_opt in inputs:
            name = phase_opt.get_name()
            output = ParsePhaseInput(phase_opt.get_dict())
            outputs += "%s: %s, " % (name, output)
    elif type(inputs) == str:
        # read file
        Utilities.my_assert(os.access(inputs, os.R_OK), FileExistsError, "Json file doesn't exist.")
        cdpt_opt = CDPT_OPT()
        cdpt_opt.read_json(inputs)
        outputs = ParsePhaseTransitionFile(cdpt_opt.toPPTF())
    else:
        # error
        raise TypeError("type of inputs must be list or str, not %s" % type(inputs))
    return outputs


def Get_reference_density(inputs):
    '''
    Get reference density of phases
    Inputs:
        inputs (dict): dictionary of phase transitions
    Returns:
        rho_p (list): reference density of phases
    '''
    # todo
    # compute density
    # read in density of the base phase
    rho0 = inputs.get("rho0", 3300.0)
    Utilities.my_assert(type(rho0) == float, TypeError, "base value for density must be a float")
    # read in the density change with each transtion
    drho = inputs.get("drho", [0.0])
    Utilities.my_assert(type(drho) == list, TypeError, "value of density change must be a list")
    total = len(drho)  # number of phase transitions
    # read in the fraction with each transtion
    xc = inputs.get("xc", [1.0])
    Utilities.my_assert(len(xc) == total, ValueError, "length of xc and drho must be the same")
    Utilities.my_assert(type(xc) == list, TypeError, "value of fraction with transition must be a list")
    rho_p = rho0 * np.ones(total + 1)
    for i in range(total):
        rho_p[i+1:] += drho[i] * xc[i]
    return rho_p
    

def Get_effective_density_change_on_660(phase_opt):
    '''
    Derive effective density change on 660 for latent heat computation
    Inputs:
        phase_opt (PHASE_OPT): options of phase transitions
    Returns:
        drho_660_eff (float): effective density change on the 660
    '''
    # todo
    lh_660 = 0.0
    name = phase_opt.get_name()
    phase_dict = phase_opt.get_dict()
    lh_660 += Get_latent_heat_contribution(phase_dict)
    rho_p = Get_reference_density(phase_dict)
    # get the index of ringwoodite phase
    i660_rd = 0
    for i in len(phase_dict['boundary']):
        if phase_dict['boundary'] == "660_ol":
            i660_rd = i + 1
            break
    rho_660 = rho_p[i660_rd] + rho_p[i660_rd+1]
    cl_660 = phase_dict['claperon slope'][i660_rd-1]
    print("rho660: %s, cl_660: %s" % (rho_660, cl_660))  # debug
    drho_660_eff = lh_660 * rho_660^2.0 / cl_660
    return drho_660_eff


def Get_latent_heat_contribution(inputs):
    '''
    Get the contribution to total latent heat from a single phase transition
    Inputs:
        inputs(dict): dictionary of a single phase transition
    Returns:
        lh_contribution (float): contribution to total latent heat
    '''
    # todo
    lh_contribution = 0.0
    return lh_contribution


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
        outputs = ParsePhaseTransitionFile(arg.inputs)
        # print the output
        print(outputs)
    elif _commend == "effective_660":
        # todo
        drho_660_eff =  0.0
        print("Effective density jump on 660 for latent heat computation is %s" % drho_660_eff)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()