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
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec 

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

class MANTLE_ADIABAT():
    '''
    Class for getting the mantle adiabatic temperature
    Attributs:
        cp - heat capacity
        alpha - thermal expansivity
        g - gravitational acceleration
        Ts = surface adiabatic temperature
    '''
    def __init__(self, cp, alpha, g, Ts, **kwargs):
        '''
        Initiation
        Inputs:
            cp - heat capacity
            alpha - thermal expansivity
            g - gravitational acceleration
            Ts = surface adiabatic temperature
            kwargs:
                approx - type of approximation
                approx_scheme - the index of this approximation
        '''
        self.cp = cp
        self.alpha = alpha
        self.g = g
        self.Ts = Ts
        self.approx = kwargs.get('approx', 'constant variables')
        if self.approx == 'constant variables':
            self.approx_scheme = 0
        elif self.approx == "constant variables and gradient":
            self.approx_scheme = 1
        else:
            raise NotImplementedError()
    
    def Temperature(self, z):
        '''
        derive temperature at depth
        Inputs:
            z - depth in m
                z could be given as value or as np.array
        '''
        if self.approx_scheme == 0:
            # cp, alpha and g are constant, approximation for the eartch condition
            T = self.Ts * np.exp(self.alpha * self.g / self.cp * z)
        elif self.approx_scheme == 1:
            # assume the thermal gradient is a constant from the surface value
            gradient = self.alpha * self.g / self.cp * self.Ts
            T = gradient * z + self.Ts
        else:
            return NotImplementedError()
        return T


class PLATE_MODEL():
    '''
    Class for the plate model temperature
    Attributes:
        L - the maximum depth of the plate model
        kappa - thermal diffusivity
        Ttop - temperature along the top boundary
        Tbot - temperature along the bottom boundary
        u - the spreading velocity
        lateral_variation - Bool, whether lateral variation
            is included in the model. If a value of u is
            given when initiating, this is set to True. Otherwise
            it is set to false.
    '''
    def __init__(self, L, kappa, Ttop, Tbot, u=None):
        '''
        Initiation
        Inputs:
            L - the maximum depth of the plate model
            kappa - thermal diffusivity
            Ttop - temperature along the top boundary
            Tbot - temperature along the bottom boundary
            u - the spreading velocity
            sommation - number of the sommation
        '''
        self.L = L
        self.kappa = kappa
        self.Ttop = Ttop
        self.Tbot = Tbot
        self.sommation = 100
        if u is None:
            self.lateral_variation = False
        else:
            self.u = u
            self.lateral_variation = True
        

    def PM_A(self, n, t):
        '''
        Get the factor A for the plate model.
        This is the factor with the \Sum term before
        the sin(n*pi*y/L) term
        Inputs:
            n - the sommation number
            t - time in s
        Returns:
            A_n(t)
        '''
        expo = (self.u * self.L / 2.0 / self.kappa -\
                (self.u**2.0 * self.L**2.0/4.0/self.kappa**2.0 + n**2.0 * np.pi**2.0)**0.5) *\
                self.u * t / self.L
        A_n_t = 2.0 / n / np.pi * \
                np.exp(expo)
        return A_n_t
    
    
    def PM_B(self, k, t):
        '''
        Get the factor B for the plate model.
        This is the factor with the \Sum term for
        the integrated heat content.
        All the even terms (n = 2k) are zero.
        Thus, here only the odd terms are handled
        Inputs:
            k - the sommation number
                note there is a relation between
                k and n:
                    n  = 2k + 1
            t - time in s
        Returns:
            B_n(t)
        '''
        n = 2*k + 1
        expo = (self.u * self.L / 2.0 / self.kappa -\
                (self.u**2.0 * self.L**2.0/4.0/self.kappa**2.0 + n**2.0 * np.pi**2.0)**0.5) *\
                self.u * t / self.L
        B_n_t = 4.0 * self.L / n**2.0 / np.pi**2.0 * \
                np.exp(expo)
        return B_n_t


    def heating_thickness(self, t):
        '''
        this is the heating thickness defined as
        heat = rho * cp * (Tp - Ts) * heating thickness
        In meaning, the heat stored in this thickness below
        the surface is released through the surface by
        conduction previous to time t in the model.
        Inputs:
            t - time to compute the heating thickness, in s
        '''
        sum = 0.0
        for i in range(self.sommation):
            sum += self.PM_B(i, t)
        heating_thickness_plate = self.L / 2.0 - sum
        return heating_thickness_plate

    def T(self, y, t):
        '''
        Get the temperature from the plate model
        todo
        '''
        if self.lateral_variation is True:
            raise NotImplementedError
        else:
            # T&S chapter 4.17
            # T0
            if type(y) == np.ndarray:
                T = np.zeros(y.size)
            else:
                T = 0.0
            T += self.Ttop
            # (T1 - T0) * y / yL0
            T += (self.Tbot - self.Ttop) * y / self.L
            for n in range(1, self.sommation):
                expo = - self.kappa * n**2.0 * np.pi**2.0 * t / self.L**2.0
                sino = n * np.pi * y / self.L
                T += (self.Tbot - self.Ttop) * 2 / np.pi / n * np.exp(expo) * np.sin(sino)
            return T



def PlotResidueHeatRatio(max_depth_plate_model, kappa, u):
    '''
    plot the ratio for the residue heat contents between the plate
    model and the half space cooling model
    Inputs:
        max_depth_plate_model - the maximum depth in the plate model (m)
        kappa - the thermal diffusivity
        u - velocity in the plate model, here in m/yr
    '''
    # todo_residue
    year = 365 * 24 * 3600.0  # s in year
    n_t = 1000 # length of the t array
    sommation = 100 # number of terms in the sum of the plate model
    Ts = 273.0
    Tp = 1673.0 # surface and mantle potential temperature
    t_array = np.linspace(0.0*year, 200e6*year, n_t)
    heating_thickness_hpc_array = np.zeros(t_array.shape)
    heating_thickness_plate_array = np.zeros(t_array.shape)
    heating_thickness_hpc_array_warm_surface = np.zeros(t_array.shape)
    heating_thickness_plate_array_warm_surface = np.zeros(t_array.shape)
    # intiating the plate model
    # the value of the surface (273.0) and the bottom (1673.0) doesn't
    # affect the calculation here
    Pmodel = PLATE_MODEL(max_depth_plate_model, kappa, Ts, Tp, u/year)
    # initiating the plot
    fig = plt.figure(tight_layout=True, figsize=(5, 10)) 
    gs = gridspec.GridSpec(2, 1)
    # plot 1: heating thickness when the surface temperature is held constant
    ax = fig.add_subplot(gs[0, 0]) 
    for i_t in range(n_t):
        t = t_array[i_t]
        heating_thickness_hpc = 2.0 * (kappa * t / np.pi)**0.5
        heating_thickness_hpc_array[i_t] = heating_thickness_hpc
        heating_thickness_plate = Pmodel.heating_thickness(t)
        heating_thickness_plate_array[i_t] = heating_thickness_plate
    ax.plot(t_array / (1e6*year), heating_thickness_hpc_array/1e3, color='tab:blue', label="hpc")
    ax.plot(t_array / (1e6*year), heating_thickness_plate_array/1e3, color='tab:red', label="plate")
    ax.set_xlabel('Age (Ma)')
    ax.set_ylabel('Heating Thickness (km)')
    # plot 2: heating thickness with a varied surface temperature
    # here, the variable t_ab stands for the age of the subducting slab
    # and t is an effective age on the slab, where is minimum temperature
    # is hotter than the surface temperature of the model as the slab
    # is warming up.
    t_ab = 40e6 * year
    t_tip = 44e6 * year
    T_minimum_at_tip = Ts + 0.2 * (Tp - Ts)
    t_array_1 = np.linspace(t_ab, t_tip, n_t)
    T_minimum_array = T_minimum_at_tip * (t_array_1 - t_ab)/ (t_tip - t_ab) + Ts * (t_array_1 - t_tip)/ (t_ab - t_tip)
    for i_t in range(n_t):
        t = t_array_1[i_t]
        Tm = T_minimum_array[i_t]
        heating_thickness_hpc_ab = 2.0 * (kappa * t_ab / np.pi)**0.5
        heating_thickness_hpc = 2.0 * (kappa * t / np.pi)**0.5
        heating_thickness_hpc_array_warm_surface[i_t] = \
            (Tm - Tp) / (Ts - Tp) * heating_thickness_hpc - heating_thickness_hpc_ab
        heating_thickness_plate_ab = Pmodel.heating_thickness(t_ab)
        heating_thickness_plate = Pmodel.heating_thickness(t)
        heating_thickness_plate_array_warm_surface[i_t] = \
            (Tm - Tp) / (Ts - Tp) * heating_thickness_plate - heating_thickness_plate_ab
    ax = fig.add_subplot(gs[1, 0]) 
    ax.plot(t_array_1 / (1e6*year), heating_thickness_hpc_array_warm_surface/1e3, color='tab:blue', label="hpc")
    ax.plot(t_array_1 / (1e6*year), heating_thickness_plate_array_warm_surface/1e3, color='tab:red', label="plate")
    ax.set_xlabel('Age (Ma)')
    ax.set_ylabel('Heating Thickness (km)')
        

    # save figure
    fig_path = os.path.join(ASPECT_LAB_DIR, "results", "residue_heat_ratio_hpc_vs_plate.png")
    fig.savefig(fig_path)
    print("%s: Generate %s" % (Utilities.func_name(), fig_path))
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
    elif _commend == 'plot_residue_heat_ratio':
        # todo_residue
        # plot the ratio for the residue heat contents between the plate
        # model and the half space cooling model
        PlotResidueHeatRatio(150e3, 1e-6, 0.05)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()