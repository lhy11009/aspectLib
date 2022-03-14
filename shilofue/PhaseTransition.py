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
        self.add_key("Clapeyron slop", list, ["clapeyron slope"], [], nick='cl')
        self.add_key("Name of boundaries", list, ["boundary"], [], nick='boundary')
        self.add_key("Depth", list, ["depth"], [], nick="depth")
        self.add_key("Width", list, ["width"], [], nick="width")
        self.add_key("Temperature", list, ["temperature"], [], nick="temperature")

    def check(self):
        '''
        check
        '''
        Utilities.JSON_OPT.check(self) # call parental function
        # check the length of vectors
        drho = self.values[1]
        xc = self.values[2]
        cl = self.values[4]
        boundary = self.values[5]
        depth = self.values[6]
        width = self.values[7]
        temperature = self.values[8]
        length = len(drho)
        assert(length == len(xc))
        assert(length == len(cl))
        if boundary == []:  # assign default values bd1 to bdn
            for i in range(length):
                boundary.append("bd%d" % (i+1))
        assert(length == len(boundary))
        assert(length == len(depth))
        assert(length == len(width))
        assert(length == len(temperature))

    def get_name(self):
        '''
        get name
        return: 
            name (str)
        '''
        return self.values[3]

    def get_boundary_names(self):
        return self.values[5]
    
    def get_depth(self):
        return self.values[6]

    def get_width(self):
        return self.values[7]
    
    def get_temperature(self):
        return self.values[8]
    
    def get_clapeyron_slope(self):
        return self.values[4]

    def to_get_reference_density(self):
        '''
        todo
        interface to get_reference_density()
        '''
        return self.values[0], self.values[1], self.values[2]

    def to_get_latent_heat_contribution(self):
        return  self.values[0], self.values[1], self.values[2], self.values[4]
    
    def total(self):
        '''
        get total number of phase transitions
        '''
        return len(self.values[1])
    
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
            "name": self.values[3],
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
    
    def get_compositions(self):
        '''
        An interface to the ParsePTfromJson
        '''
        return self.values[0]


def ParsePhaseInput(phase_opt, entry="density"):
    '''
    parse input of phases to a aspect input form
    Inputs:
        entry (str): one in ["density", ], indicating what type of output we want
    
    '''
    output = ""
    # number of transition
    total = phase_opt.total()
    if entry == "density":
        rho = Get_reference_density(*phase_opt.to_get_reference_density())
        # generate output
        output += "%.1f" % rho[0]
        for i in range(1, total+1):
            output += "|%.1f" % rho[i]  
    elif entry == "depth":
        depth = phase_opt.get_depth()
        # generate output
        output += "%.1f" % depth[0]
        for i in range(1, total):
            output += "|%.1f" % depth[i]  
        pass
    elif entry == "width":
        width = phase_opt.get_width()
        # generate output
        output += "%.1f" % width[0]
        for i in range(1, total):
            output += "|%.1f" % width[i]  
        pass
    elif entry == "temperature":
        temperature = phase_opt.get_temperature()
        # generate output
        output += "%.1f" % temperature[0]
        for i in range(1, total):
            output += "|%.1f" % temperature[i]  
        pass
    elif entry == "clapeyron slope":
        cl = phase_opt.get_clapeyron_slope()
        # generate output
        output += "%.1f" % cl[0]
        for i in range(1, total):
            output += "|%.1f" % cl[i]  
        pass
    return output


def ParsePTfromJson(inputs):
    '''
    Parse the inputs of phase transition from a json file to output in aspect
    Inputs:
        inputs (list): a list of phase transitions (phase_opt objects)
                or
                (str): path to a json file
    Returns:
        outputs (dict): a dictionary for a aspect prm file
    '''
    outputs = "" 
    if type(inputs) == list:
        # read dict
        # get the outputs
        outputs = {}
        density_outputs = ""
        # density output
        is_first = True
        for phase_opt in inputs:
            name = phase_opt.get_name()
            # todo
            output = ParsePhaseInput(phase_opt)
            if is_first:
                is_first = False
            else:
                density_outputs += ', '
            density_outputs += "%s: %s" % (name, output)
        outputs["Densities"] = density_outputs
        # depth
        is_first = True
        depth_outputs = ""
        for phase_opt in inputs:
            name = phase_opt.get_name()
            output = ParsePhaseInput(phase_opt, "depth")
            if is_first:
                is_first = False
            else:
                depth_outputs += ', '
            depth_outputs += "%s: %s" % (name, output)
        outputs["Phase transition depths"] = depth_outputs 
        # width    
        is_first = True 
        width_outputs = ""
        for phase_opt in inputs:
            name = phase_opt.get_name()
            output = ParsePhaseInput(phase_opt, 'width')
            if is_first:
                is_first = False
            else:
                width_outputs += ', '
            width_outputs += "%s: %s" % (name, output)
        outputs["Phase transition widths"] = width_outputs 
        # temperature
        is_first = True
        temperature_outputs = ""
        for phase_opt in inputs:
            name = phase_opt.get_name()
            output = ParsePhaseInput(phase_opt, 'temperature')
            if is_first:
                is_first = False
            else:
                temperature_outputs += ', '
            temperature_outputs += "%s: %s" % (name, output)
        outputs["Phase transition temperatures"] = temperature_outputs
        # clapeyron slope
        is_first = True
        cl_outputs = ""
        for phase_opt in inputs:
            name = phase_opt.get_name()
            output = ParsePhaseInput(phase_opt, 'clapeyron slope')
            if is_first:
                is_first = False
            else:
                cl_outputs += ', '
            cl_outputs += "%s: %s" % (name, output)
        outputs["Phase transition Clapeyron slopes"] = cl_outputs
    elif type(inputs) == str:
        # read file
        Utilities.my_assert(os.access(inputs, os.R_OK), FileExistsError, "Json file doesn't exist.")
        cdpt_opt = CDPT_OPT()
        cdpt_opt.read_json(inputs)
        outputs = ParsePTfromJson(cdpt_opt.get_compositions())
    else:
        # error
        raise TypeError("type of inputs must be list or str, not %s" % type(inputs))
    return outputs


def Get_reference_density(rho0, drho, xc):
    '''
    Get reference density of phases
    Inputs:
    Returns:
        rho_p (list): reference density of phases
    '''
    # todo
    # compute density
    Utilities.my_assert(type(rho0) == float, TypeError, "base value for density must be a float")
    Utilities.my_assert(type(drho) == list, TypeError, "value of density change must be a list")
    total = len(drho)  # number of phase transitions
    # read in the fraction with each transtion
    Utilities.my_assert(len(xc) == total, ValueError, "length of xc and drho must be the same")
    Utilities.my_assert(type(xc) == list, TypeError, "value of fraction with transition must be a list")
    rho_p = rho0 * np.ones(total + 1)
    for i in range(total):
        rho_p[i+1:] += drho[i] * xc[i]
    return rho_p


def Get_entropy_change(rho0, drho, xc, cl):
    '''
    Get entropy changes on individual phases
    Equation used here:
        dS = - gamma * drho / rho^2.0
    Returns:
        lh_contribution (float): contribution to total latent heat
    '''
    # todo
    lh=[]
    rho_p = Get_reference_density(rho0, drho, xc)
    for i_dx in range(len(drho)):
        rho = (rho_p[i_dx] + rho_p[i_dx+1]) / 2.0
        lh.append( -1.0 * cl[i_dx] * drho[i_dx] *xc[i_dx] / rho**2.0)
    return lh


def Get_temperature_change(rho0, drho, xc, cl, **kwargs):
    '''
    Get entropy changes on individual phases
    Equation used here:
        dS = - gamma * drho / rho^2.0
        dT / T0 = 1 / (1 + Ds / cp)
    kwargs (dict):
        cp - heat capacity
    Returns:
        dT: factors of temperature variation from latent heat
    '''
    # todo
    cp = kwargs.get("cp", 1250.0)
    dS = Get_entropy_change(rho0, drho, xc, cl)
    dT = []
    for i in range(len(dS)):
        dT.append(1.0/(1.0+dS[i]/cp))
    return dT
        


def Show_entropy_changes(file_path):
    '''
    todo
    Use the Get_entropy_change function and show entropy changes on phase transitions
    Inputs:
        file_path(str): a configuration for phases
    '''
    assert(os.access(file_path, os.R_OK))
    cdpt_opt = CDPT_OPT()
    cdpt_opt.read_json(file_path)
    compositions = cdpt_opt.get_compositions()
    for i in range(len(compositions)):
        phase_opt = compositions[i]
        dS = Get_entropy_change(*phase_opt.to_get_latent_heat_contribution())
        print("Composition: ", phase_opt.get_name())
        print("Entropy changes [J/(Kg*K)]")
        print(phase_opt.get_boundary_names())  # screen output
        print(dS)
        dT = Get_temperature_change(*phase_opt.to_get_latent_heat_contribution())
        print("Factors of temperature change (T2/T0)")
        print(phase_opt.get_boundary_names())  # screen output
        print(dT)
        print('\n')


def Usage():
    CDPT_opt=CDPT_OPT()
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
    - parse composition dependent phase transition inputs:\n\
\n\
        python -m shilofue.PhaseTransition phase_input -i ./files/TwoDSubduction/phases_1_0.json\n\
\n\
    - show the entropy changes on phase transitions\n\
\n\
        python -m shilofue.PhaseTransition show_entropy_changes -i /home/lochy/ASPECT_PROJECT/aspectLib/tests/integration/fixtures/phase_transitions/phases.json\n\
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
        outputs = ParsePTfromJson(arg.inputs)
        # print the output
        print(outputs)
    elif _commend == "show_entropy_changes":
        Show_entropy_changes(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()