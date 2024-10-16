# -*- coding: utf-8 -*-
r"""hk03 rheology and ploting, I got this from Magali Billen in 2021

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
    Originally, it only has diffusion creep rheology
""" 

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf


ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')

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
  - plot a PT map for the peierls rheology: \n\
\n\
        python -m shilofue.FlowLaws plot_peierls_creep_PT -f MK10\n\
\n\
        -f: type of flow law to use (MK10 or Idrissi16, the second one doesn't converge for now)\n\
\n\
  - plot a deformation mechanism for peierls rheology: \n\
\n\
        python -m shilofue.FlowLaws plot_peierls_dfm\n\
\n\
        ")


#Physical Constants
R=8.314 #J/mol*K


# Input parameter choices as:
# T = temperature               (K)
# P = pressure                (Pa)
# d = grainsize               (microns) -> converted to meters inside.
# sigd = differential stress  (MPa)
# coh = water content         (ppm)
# fh2o = MPa

# Edev, Vdev:  "mid", "min" or "max" value from reported error range.
# mod:  "orig" or "new" refers to using the orignal values of V (and A) from HK03 or new values 
# following Ohuchi 2015.

# Most flow laws require the water fugacity rather than the water content, but we usually
# define water as C_OH, so need to convert this to pass to flow-laws.
# ** Exception is the original HK03 "Constant C_OH" flow laws which already did this conversion.
 
def convert_OH_fugacity(T,P,coh):

    # T [K]             Temperature
    # P [MPa]            Pressure
    # coh [ppm-H/si]    OH content (note can convert C_h2o = Coh/16.25, confirmed by Ohuchi)
    
    # Per e-mail from Ohuchi, use data from the table in Keppler & Bolfan-Casanova, 2006
    # to related C_H2O to f_H2O, checked against Kohlstedt et al 1996 data. 
    # I used Kohlstedt et al., 1996 data to show that CH2O = COH/16.25
    # substituting this leads to multiply Ah2o by 16.25Then change units
    # bars to MPa... the resulting value (87.75 ppm-H/Si/MPAa) is very close to the 
    # Ah2o in Zhao (90 # ppm-H/Si/MPa +/- 10) but we ignore the iron content dependence.  
       
    Ah2o = 0.54*(10)*16.25             # (ppm/bar)*(MPa/bar)*COH/CH2O
    Eh2o = 50e3                     # J/mol +/-2e3
    Vh2o = 10.6e-6                     # m^3/mol+/-1
    fH2O = coh/(Ah2o*np.exp(-(Eh2o + P*Vh2o)/(R*T)))       #MPa
    
    return fH2O

# Use choices of water (wet/dry/constant), mod (orig,modified) and deviation (mid, min, max)
# to set the values of A, E and V to use in the strain rate or viscosity calculation
# functions.

# ** I think this information might be better kept in a dictionary?? ** 
# wet or dry or con
# uses "mid" reported values, 
# and define dE or dV, which is +/- below
       
def get_diff_HK_params(water,mod,Edev,Vdev):
    p = 3 
    
    if water == 'dry':
        r = 0
        if mod == 'orig':
            A = 1.5e9/(1e6)**p # MPa^(-n-r) * m^p * s^-1 (converted from microns to meters)
            E = 375e3      # J/mol 
            V = 6e-6       # m^3/mol, Mid value form range in Table 1 of HK 03.
            dE = 50e3      # error on activation energy
            dV = 4e-6      # error on activation volume 
        else: 
            A = 10**(-10.40) #MPa^(-n-r) * m^p * s^-1  Hansen et al (2011), Table A2 logAdif=7.6 MPa, microns, K; grain size difference; convert microns to meters 
            E = 375e3       # J/mol     
            V = 8.2e-6      # m^3/mol       Nishihara et al. (2014) for forsterite   8.2 +/- 0.9 cm^3.mol
            dE = 50e3          # error on activation energy from HK03
            dV = 0.9e-6      # error on activation volume from Nishihara         
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1
        if mod == 'orig':
            A = 10**(-10.6) # MPa^(-n-r) * m^p * s^-1 # HK03: 2.5e7 converted from microns to meters        
            E = 375e3         # J/mol                                    
            V = 10e-6          # m^/mol, from HK03   
            dE = 75          # error on activation energy, from HK03
            dV = 10e-6      # error on activation volume, from HK03  
        else: 
            A = (10**(-10.6))/3.5 #MPa^(-n-r) * m^p * s^-1 # 2.5e7 converted from microns to meters; bell03 water correction        
            E = 375e3  # J/mol                                    
            V = 23e-6  # m^/mol, Ohuchi et al. (2012)   
            dE = 75e3    # error on activation energy from HK03
            dV = 4.5e-6  # Ohuchi et al, 2012
    elif water == 'con': 
        print('Using Diffusion Constant-COH')     
        r = 1        
        A = 1.0e6/(1e6)**p # Pa^(-n-r) * m^p * s^-1 (converted from microns to meters)      
        E = 335e3         # J/mol                                    
        V = 4e-6          # m^/mol, from HK03   
        dE = 75e3          # error on activation energy, from HK03
        dV = 4e-6          # use error for from HK03  
        if mod != 'orig':
            print("Error no modified versions for Constant-COH HK03, using originals")
    else:
        print("water must be dry, wet or con")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
        
    return(A,E,V,p,r)
    
    
def get_disl_HK_params(water,mod,Edev,Vdev):
    p = 0   # no grain size dependence
    
    if water == 'dry':
        r = 0
        n = 3.5
        dn = 0.3
        if mod == 'orig':
            A = 1.1e5 # MPa^(-n-r) * s^-1 (=10^5.04)
            E = 530e3      # J/mol 
            V = 18e-6   # m^3/mol, Mid value from range in Table 2 of HK 03.
            dE = 4e3      # error on activation energy
            dV = 4.0e-6      # error on activation volume 
        else: 
            A = 1.1e5 # MPa^(-n-r) * s^-1 (=10^5.04)
            E = 530e3      # J/mol 
            V = 17e-6   # m^3/mol, from Kawazoe etal., 2009            
            dE = 4e3      # error on activation energy
            dV = 2.5e-6      # error on activation volume     
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1.2
        dr = 0.4
        n = 3.5
        dn = 0.3
        if mod == 'orig':
            A = 1600         # MPa^(-n-r) # HK03 (=10^3.2)
            E = 520e3         # J/mol                                    
            V = 22e-6          # m^/mol, from HK03   
            dE = 40e3          # error on activation energy, from HK03
            dV = 11e-6      # error on activation volume, from HK03  
        else: 
            A = 1600/3.5 #MPa^(-n-r) s^-1 # HK03 modified for Bell_etal 03: (=10^2.66)
            E = 520e3  # J/mol                                    
            V = 24e-6  # m^/mol, Ohuchi et al. (2012)   
            dE = 40e3    # error on activation energy from HK03
            dV = 3e-6  # NEED TO CHECK Ohuchi et al, 2012
    elif water == 'con': 
        print('Using Dislocation Constant-COH')     
        r = 1.2 
        dr = 0.4
        n = 3.5
        dn = 0.3       
        A = 90 # MPa^(-n-r) * s^-1     
        E = 480e3         # J/mol                                    
        V = 11e-6          # m^/mol, from HK03   
        dE = 40e3          # error on activation energy, from HK03
        dV = 5e-6          # use error for from HK03  
        if mod != 'orig':
            print("Error no modified versions for Constant-COH HK03, using originals")
    else:
        print("water must be dry, wet or con")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
        
    return(A,E,V,n,r)

# Preferred dislocation creep parameters from Kawazoe 
def get_disl_KAW09_params(water,mod,Edev,Vdev):
    p = 0   # no grain size dependence
    
    if water == 'dry':  # same as modified HK03
        r = 0
        n = 3.5
        dn = 0.3
        A = 1.1e5 # MPa^(-n-r) * s^-1 (=10^5.04)
        E = 530e3      # J/mol 
        V = 17e-6   # m^3/mol, from Kawazoe etal., 2009            
        dE = 4e3      # error on activation energy
        dV = 2.5e-6      # error on activation volume    
    elif water == 'wet':      # HK03 wet diffusion, (NOT constant COH)
        r = 1.2
        dr = 0.4
        n = 3.0
        dn = 0.1
        A = 10**2.9 #MPa^(-n-r) s^-1 # HK03 modified for Bell_etal 03: (=10^2.66)
        E = 470e3  # J/mol                                    
        V = 24e-6  # m^/mol,  
        dE = 40e3  #
        dV = 3e-6  # 
    else:
        print("water must be dry or wet")
                    
    if Edev == 'min':  
        E = E - dE  
    elif Edev == 'max':
        E = E + dE  

    if Vdev == 'min':  
        V = V - dV  
    elif Vdev == 'max':
        V = V + dV
    return(A,E,V,n,r)    
    
# Olivine diffusion creep: strain-rate
def edot_diff_HK(T,P,d,sigd,coh,water,mod,Edev,Vdev):   

    A, E, V, p, r = get_diff_HK_params(water,mod,Edev,Vdev)    

    if water == 'con': # use constant COH equation/values from HK03
        edot = A*(sigd)**n*(d/1e6)**(-p)*((coh)**r)*np.exp(-(E+P*V)/(R*T)) #s^-1
    else:
        fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
        edot = A*(sigd)**n*(d/1e6)**(-p)*((fh2o)**r)*np.exp(-(E+P*V)/(R*T)) #s^-1

    return edot
    
# Olivine dislocation creep: strain-rate
def edot_disc_HK(T,P,sigd,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disc_HK_params(water,mod,Edev,Vdev)    

    if water == 'con': # use constant COH equation/values from HK03
        edot = A*(sigd)**n*((coh)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1
    else:
        fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
        edot = A*(sigd)**n*((fh2o)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1

    return edot
    
# Olivine dislocation creep: strain-rate
def edot_disc_KAW09(T,P,sigd,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disc_KAW09_params(water,mod,Edev,Vdev)    

    fh2o = convert_OH_fugacity(T,P,coh) # Convert coh to water fugacity
    edot = A*(sigd)**n*((fh2o)**r)*np.exp(-(E+P*V)/(n*R*T)) #s^-1

    return edot

# Olivine diffusion creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
def visc_diff_HK(T,P,d,coh,water,mod,Edev,Vdev):   

    A, E, V, p, r = get_diff_HK_params(water,mod,Edev,Vdev)
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
        
    if water == 'con': # use constant COH equation/values from HK03        
        visc = 0.5*(dm)**p/(Am*coh**r)*np.exp((E + P*V)/(R*T)) # Pa s
    else:     
        fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
        visc = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((E + P*V)/(R*T)) # Pa s

    return visc
    
# Olivine dislocation creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
def visc_disl_HK(T,P,edot,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disl_HK_params(water,mod,Edev,Vdev)    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
        
    if water == 'con': # use constant COH equation/values from HK03        
        visc = 0.5*edot**(1/n - 1)/(Am*coh**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s
    else:     
        fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
        visc =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s

    return visc

# Olivine dislocation creep: viscosity  
# Note that A in this case needs to be convert from MPa to Pa to give visc in Pa-s
# instead of MPa-s 
# For Kawazoe et al., 2009 preferred values    
def visc_disl_KAW09(T,P,edot,coh,water,mod,Edev,Vdev):   
    
    A, E, V, n, r = get_disl_KAW09_params(water,mod,Edev,Vdev)    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
        
    fh2o = convert_OH_fugacity(T,P,coh)    # Convert coh to water fugacity
    visc =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((E + P*V)/(n*R*T)) # Pa s

    return visc

def plot_upper_mantle_viscosity(**kwargs):
    '''
    plot upper mantle viscosity, I got this from Magali
    '''
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = "10"
    
    # Unit conversions
    mpa = 1e6  # MPa to Pa
    mum = 1e6  # microns/meter for grain size 
    R = 8.314       # J/mol*K
    km2m = 1e3 # km to meters
    sec2yrs = 60*60*24*365.25 # sec per year
    
    # Depth in meters
    zmin = 0
    zmax = 660*km2m
    dz = 20*km2m
    depth = np.arange(zmin,zmax+dz,dz)
    
    # Temperature for Half-space Cooling Model
    age = 80e6*sec2yrs
    T = temperature_halfspace(depth,age,0) # K
    
    # temperature for flow law
    mad = 0.5715  # deg/km, approximating adiabat to match T(CMB) in Aspect
    Tad = T + mad*depth/km2m
    
    P = pressure_from_lithostatic(depth,Tad,0)  # Pa
    
    ptplt = kwargs.get('ptplt', 0)
    if ptplt == 1: 
        fig = plt.figure()
        ax = fig.add_subplot(2,3,1)
        ax.plot(P/1e9,depth/km2m,color='blue')
        ax.set_ylim(zmax/km2m,0)
        ax.grid(True)
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Depth (km)')
    
        ax = fig.add_subplot(2,3,2)
        ax.plot(T,depth/km2m,color='green')  # plot in kilometers
        ax.plot(Tad,depth/km2m,color='blue')  # plot in kilometers
        ax.set_ylim(zmax/km2m,0)
        ax.grid(True)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Depth (km)')
        
        plt.tight_layout()
        pdffile = os.path.join(RESULT_DIR, 'profile_pt.pdf')
        fig.savefig(pdffile,bbox_inches='tight')
        
    # Use for testing at a single value
    #T = 1400 + 273
    #P = 1e9
    #print(f'etad {etadf:0.3g}')
    
    # constant grain size and strain rate
    d = 5e3       # microns (= 5 mm = 0.5 cm)
    edot = 1e-15  # background mantle strain-RateLimiter
    
    # Which flow law?
    mech = 'Diffusion'
    
    # Calculate diffusion creep viscosity (water, modified?, Edev, Vdev)
    # Which flow law version?
    flver = kwargs.get('flavor', 'orig')
    
    # Use top or bottom of error range (deviation from mid?) min/mid/max
    Edeva = ['min','mid','max']
    Vdeva = ['min','mid','max']
    cnt = 1
    fig1 = plt.figure(figsize=[7.5,9.0])
    for i in range(3):
        print('Edev',Edeva[i])
        Edev = Edeva[i]
        
        for j in range(3):
            print('Vdev',Vdeva[j])    
            Vdev = Vdeva[j]
            
            titext1 = 'E: ' + Edev + ', V:' + Vdev
            # Water: choose wet or dry and indicate water content for wet
            water = 'wet'
            coh = 1000  # ppm H/Si
            etadf_wet = visc_diff_HK(Tad,P,d,coh,water,flver,Edev,Vdev)
    
            water = 'dry'
            coh = 50  # ppm H/Si
            etadf_dry = visc_diff_HK(Tad,P,d,coh,water,flver,Edev,Vdev)
    
            water = 'con'
            coh = 1000  # ppm H/Si
            etadf_con = visc_diff_HK(Tad,P,d,coh,water,flver,Edev,Vdev)
            print('temperature: %4e, pressure: %4e' % (Tad[-1], P[-1]))  # debug
            print('viscosity: %.4e' % etadf_con[-1])

            ax = fig1.add_subplot(3,3,cnt)
            ax.semilogx(etadf_wet,depth/km2m,label='wet')
            ax.semilogx(etadf_dry,depth/km2m,label='dry')
            ax.semilogx(etadf_con,depth/km2m,label='const')
            ax.set_ylim(zmax/km2m,0)
            ax.set_xlim(1e18,1e25)
            ax.grid(True)
            ax.minorticks_on()    
            ax.set_xticks([1e18, 1e19, 1e20, 1e21, 1e22, 1e23, 1e24])
            
            ax.legend(fontsize=8)
            if i == 2:
                ax.set_xlabel('Viscosity (Pa s)')
            if j == 0:
                ax.set_ylabel('Depth (km)')
            if cnt == 1:
                titext2 = mech + '(' + flver + ') ' + titext1
                ax.set_title(titext2, fontsize=10)
            else:
                ax.set_title(titext1, fontsize=10)
            cnt = cnt + 1
            
    
    plt.tight_layout()
    pdffile = os.path.join(RESULT_DIR, mech + '_' + flver + '.pdf')
    fig1.savefig(pdffile,bbox_inches='tight')


# Half-space cooling model 
# t in seconds
# z in meters
def temperature_halfspace(z,t,p):
    # Physical constants
    kappa = 1e-6  # thermal diffusivity (m^2/s)
    T_s = 273  # surface temperature (K)
    T_m = 1673 # mantle temperature (K)

    T = T_s + (T_m - T_s)*erf(z/(2*np.sqrt(kappa*t)))
    
    if p == 1:
        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        ax.plot(T,z/1000)  # plot in kilometers
        ax.set_ylim(150,0)
        ax.grid(True)
        ax.set_xlabel('Temperature (C)')
        ax.set_ylabel('Depth (km)')
    
    return T


# Pressure as a function of depth from the compressibility using Equation 4-321 
# (Turcotte and Schubert, Geodynamics 2nd edition, p. 190). 
# z in meters
def pressure_from_compessibility(z,p):
    # Physical constants
    beta = 4.3e-12  # compressiblity (1/Pa)
    rho_m = 3300    # reference density (kg/m^3)
    g = 9.81        # gravitational acceleration (m/s^2)
    
    P = (-1/beta)*np.log(1 - rho_m*g*beta*z)  # pressure (Pa)
    if p == 1:
        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        ax.plot(P/1e9,z/1000,color='blue')
        ax.set_ylim(150,0)
        ax.grid(True)
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Depth (km)')
    return P


# Pressure profile from lithostatic density
# Aspect uses total pressure, so this includes a lithostatic pressure
# even if using incompressible formulation
# This requires the temperature as a function of depth (including adiabatic gradient)

def pressure_from_lithostatic(z,Tad,p):
    
    # Density Profile
    refrho = 3300  # kg/m^3
    refT = 1673        # K
    alpha = 3.1e-5  # 1/K
    g = 9.81 # m/s^2
    density = refrho*(1-alpha*(Tad-refT))

    # start loop at 1 because P[0] = 0
    dz = z[1]-z[0]
    P = np.zeros(np.size(z))
    for i in range(1, np.size(z)):
        P[i] = P[i-1] + 0.5*(density[i]+density[i-1])*g*dz
    
    if p == 1:
        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        ax.plot(P/1e9,z/1000,color='blue')
        ax.set_ylim(150,0)
        ax.grid(True)
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Depth (km)')
    return P


### peierls creep

def GetPeierlsApproxVist(flv):
    '''
    export Peierls rheology for approximation
    '''
    mpa = 1e6  # MPa to Pa
    Peierls={}
    if flv == "MK10":
        Peierls['q'] = 1.0
        Peierls['p'] = 0.5
        n = 2.0
        Peierls['n'] =  n
        Peierls['sigp0'] = 5.9e9					# Pa (+/- 0.2e9 Pa)
        Peierls['A'] = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
        Peierls['E'] = 320e3  					# J/mol (+/-50e3 J/mol)
    elif flv == "Idrissi16":
        Peierls['q'] = 2.0
        Peierls['p'] = 0.5
        n = 0.0
        Peierls['n'] =  n
        Peierls['sigp0'] = 3.8e9					# Pa (+/- 0.7e9 Pa)
        Peierls['A'] = 1e6 	# s^-1, note unit of A is related to n (here n = 0.0)
        Peierls['E'] = 566e3  					# J/mol (+/-74e3 J/mol)
    else:
        raise ValueError("flv must by \'MK10\'")
    return Peierls 


def GetPeierlsStressPDependence(flv):
    '''
    export P dependence for the peierls creep
    '''
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    return G0, Gp


def peierls_visc_from_stress(flv, P, T, sigma, **kwargs):
    '''
    Peierls creep flow law 
    flv: flow law version
    '''
    mpa = 1e6  # MPa to Pa
    if flv == "MK10":
        # Mei et al., JGR 2010
        q = 1.0
        p = 0.5
        n = 2.0
        sigp0 = 5.9e9					# Pa (+/- 0.2e9 Pa)
        A = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
        E = 320e3  					# J/mol (+/-50e3 J/mol)
        V = 0.0
    elif flv == "Idrissi16":
        # Idrissi et al 2016, this doesn't actually converge.
        q = 2.0
        p = 0.5
        n = 0
        sigp0 = 3.8e9
        A = 1e6
        E = 566e3
        V = 0.0
    else:
        raise ValueError("flv must be \"MK10\" or \"Idrissi16\".")
    prescribe_A = kwargs.get("A", None)
    if prescribe_A is not None:
        A = prescribe_A
    expo = np.exp(-(E + P*V) / (R*T) * (1 - (sigma/sigp0)**p)**q)
    edot = A * expo * sigma ** n
    eta = sigma / (2 * edot)
    return eta, edot


def peierls_visc_from_edot(flv, P, T, edotp0, limit=0.1):
    '''
    Get value of peierls creep from strain rate
    flv: flow law version
        MK10: for Mei and Kohlstedt, 2010     (gam = 0.17)
    P: pressure
    T: temperature
    edotp0: strain rate
    limit: limit of error
    '''
    diff = 1e6  # a big initial value
    sigma_l = 1e-5
    sigma_u = 1e12
    is_first = True
    n = 0
    while (abs(diff) > limit):
        if is_first:
            is_first = False
        else:
            # update the value of sigma
            if diff > 0.0:
                sigma_u = sigma
            else:
                sigma_l = sigma
        exponential = (np.log(sigma_u) + np.log(sigma_l)) / 2.0
        sigma = np.exp(exponential)
        etap, edotp = peierls_visc_from_stress(flv, P, T, sigma)
        diff = np.log(edotp / edotp0)
        n += 1
    return etap, sigma, diff, n


def peierls_visc_from_edot_newton(flv, P, T, edotp0, limit=0.1, **kwargs):
    '''
    Get value of peierls creep from strain rate, using the newton method
    flv: flow law version
        MK10: for Mei and Kohlstedt, 2010     (gam = 0.17)
    P: pressure
    T: temperature
    edotp0: strain rate
    limit: limit of error
    '''
    debug = kwargs.get("debug", False)
    i_max = kwargs.get('i_max', 20)
    mpa = 1e6  # MPa to Pa
    prescribe_A = kwargs.get("A", None)
    if flv == "MK10":
        # Mei et al., JGR 2010
        q = 1.0
        p = 0.5
        n = 2.0
        sigp0 = 5.9e9					# Pa (+/- 0.2e9 Pa)
        A = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
        E = 320e3  					# J/mol (+/-50e3 J/mol)
        gam = 0.17    # gam: fitting parameter = sig_ref/sigp
    elif flv == "Idrissi16":
        q = 2.0
        p = 0.5
        n = 0.0
        sigp0 = 3.8e9					# Pa (+/- 0.7e9 Pa)
        A = 1e6 	# s^-1, note unit of A is related to n (here n = 0.0)
        E = 566e3  					# J/mol (+/-74e3 J/mol)
        gam = 0.15  # use the same value they used in the test
    else:
        raise ValueError("type of flow law has to be \"MK10\" or \"Idrissi16\"")


    diff = 1e6  # a big initial value
    sigma_l = 1e-5
    sigma_u = 1e12
    is_first = True
    i = 0
    # using the approx value as initial guess
    visc0 = peierls_approx_visc(flv, P,T,edotp0, debug=debug, A=prescribe_A)
    if debug:
        print("visc0 = ", visc0)
    sigma = 2 * visc0 * edotp0
    _, edotp = peierls_visc_from_stress(flv, P, T, sigma, A=prescribe_A)
    diff = np.log(edotp / edotp0)
    if debug:
        print("\ni = %d, sigma = %.4e, edotp = %.4e, diff = %.4e (log (edotp / edotp0))"\
                % (i, sigma, edotp, diff))
    while (abs(diff) > limit):
        if i > i_max:
            print("peierls_visc_from_edot_newton: maximum iteration (%d) is reached" % i_max)
            break
        # compute the gradient of ln (strain_rate) - ln (sigma)
        ln_grad = n + p * q * (E / R / T) *\
                  (sigma / sigp0)**p * (1 - (sigma / sigp0)**p)**(q-1)
        # compute the new value of stress using the newton method 
        sigma = np.exp(np.log(sigma) - diff / ln_grad)
        # compute the new viscosity and strain rate
        etap, edotp = peierls_visc_from_stress(flv, P, T, sigma, A=prescribe_A)
        diff = np.log(edotp / edotp0)
        i += 1
        if debug:
            print("i = %d, sigma = %.4e, edotp = %.4e, ln_grad = %.4e, diff = %.4e (log (edotp / edotp0))"\
                     % (i, sigma, edotp, ln_grad, diff))
    return etap, sigma, diff, i


def peierls_visc_from_edot_newton_nolog(flv, P, T, edotp0, limit=0.1, **kwargs):
    '''
    Get value of peierls creep from strain rate, using the newton method, 
    compute the derivative without the log value
    flv: flow law version
        MK10: for Mei and Kohlstedt, 2010     (gam = 0.17)
    P: pressure
    T: temperature
    edotp0: strain rate
    limit: limit of error
    kwargs(dict):
        debug: debug options
    '''
    debug = kwargs.get("debug", False)
    mpa = 1e6  # MPa to Pa
    if flv == "MK10":
        # Mei et al., JGR 2010
        q = 1.0
        p = 0.5
        n = 2.0
        sigp0 = 5.9e9					# Pa (+/- 0.2e9 Pa)
        A = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
        E = 320e3  					# J/mol (+/-50e3 J/mol)
        gam = 0.17    # gam: fitting parameter = sig_ref/sigp
    elif flv == "Idrissi16":
        q = 2.0
        p = 0.5
        n = 0.0
        sigp0 = 3.8e9					# Pa (+/- 0.7e9 Pa)
        A = 1e6 	# s^-1, note unit of A is related to n (here n = 0.0)
        E = 566e3  					# J/mol (+/-74e3 J/mol)
        gam = 0.15  # use the same value they used in the test
    else:
        raise ValueError("type of flow law has to be \"MK10\" or \"Idrissi16\"")

    diff = 1e6  # a big initial value
    sigma_l = 1e-5
    sigma_u = 1e12
    is_first = True
    i = 0
    # using the approx value as initial guess
    visc0 = peierls_approx_visc(flv, P,T,edotp0)
    sigma = 2 * visc0 * edotp0
    _, edotp = peierls_visc_from_stress(flv, P, T, sigma)
    diff = np.log(edotp / edotp0)
    diff_1 = edotp - edotp0
    if debug:
        print("i = %d, sigma = %.4e, edotp = %.4e, diff = %.4e (log (edotp / edotp0))"\
                % (i, sigma, edotp, diff))
    while (abs(diff) > limit):
        # compute the gradient of ln (strain_rate) - ln (sigma)
        grad = A * sigma**(n-1) * np.exp((- E/R/T) * (1 - (sigma / sigp0)**p)**q) *\
                ( (E/R/T) * p * q * (sigma / sigp0)**p * (1 - (sigma / sigp0)**p)**(q-1) + n)
        # compute the new value of stress using the newton method
        sigma = sigma - diff_1 / grad
        # compute the new viscosity and strain rate
        etap, edotp = peierls_visc_from_stress(flv, P, T, sigma)
        diff = np.log(edotp / edotp0)
        diff_1 = edotp - edotp0
        i += 1
        if debug:
            print("i = %d, grad = %.4e, sigma = %.4e, edotp = %.4e, diff = %.4e (log (edotp / edotp0))"\
                    % (i, grad, sigma, edotp, diff))
    return etap, sigma, diff, i


def peierls_approx_visc(flv, P,T,edot, **kwargs):
    '''
    Peierls creep flow law 
    flv: flow law version
    MK10: for Mei and Kohlstedt, 2010     (gam = 0.17)
    '''
    debug = kwargs.get("debug", False)
    mpa = 1e6  # MPa to Pa

    if flv == "MK10":
        # Mei et al., JGR 2010
        q = 1.0
        p = 0.5
        n = 2.0
        sigp0 = 5.9e9					# Pa (+/- 0.2e9 Pa)
        A = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
        E = 320e3  					# J/mol (+/-50e3 J/mol)
        gam = 0.17    # gam: fitting parameter = sig_ref/sigp
    elif flv == "Idrissi16":
        q = 2.0
        p = 0.5
        n = 0.0
        sigp0 = 3.8e9					# Pa (+/- 0.7e9 Pa)
        A = 1e6 	# s^-1, note unit of A is related to n (here n = 0.0)
        E = 566e3  					# J/mol (+/-74e3 J/mol)
        gam = 0.15  # use the same value they used in the test
    else:
        raise ValueError("type of flow law has to be \"MK10\" or \"Idrissi16\"")

    # prescribe A
    prescribe_A = kwargs.get('A', None)
    if prescribe_A is not None:
        A = prescribe_A
    
    # Pressure dependence of Peierls stress
    # From Kawazoe et al. PEPI 2009 and parameters from Liu et al., GRL 2005 
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    sigp = sigp0*(1 + (Gp/G0)*P)
    
    s = (E/(R*T))*p*q*((1-gam**p)**(q-1))*(gam**p)
    x = 1/(s+n)
    visc = (0.5*gam*sigp*edot**(x-1))/( ((A*(gam*sigp)**n)**x)*np.exp( -(E*(1-gam**p)**q)/(R*T*(s+n)) ) )
    stress_term = (0.5*gam*sigp)/( ((A*(gam*sigp)**n)**x))
    strain_rate_term = edot**(x-1)
    arrehius_term = np.exp((E*(1-gam**p)**q)/(R*T*(s+n)) )
    if debug:
        print("stress_term = ", stress_term)
        print("strain_rate_term = ", strain_rate_term)
        print("arrehius_term = ", arrehius_term)
    
    return visc


### plot utilities

def PlotPeierlsPTMap(flv, implementation, edotp0):
    '''
    plot the peierls creep in P and T space
    '''
    assert(implementation in ['exact', 'approximation'])  # two types of implementation to use
    nP = 100
    nT = 200
    Ps = np.linspace(0, 2.6e10, nP)
    Ts = np.linspace(273.0, 2000, nT)
    PPs, TTs = np.meshgrid(Ps, Ts)
    etapps = np.zeros(PPs.shape)
    for i in range(PPs.shape[0]):
        for j in range(PPs.shape[1]):
            P = PPs[i, j]
            T = TTs[i, j]
            if implementation == 'exact':
                etap, _, _, _ = peierls_visc_from_edot(flv, P, T, edotp0, limit=0.1)
            elif implementation == 'approximation':
                etap = peierls_approx_visc(flv, P,T,edotp0)
            etapps[i, j] = etap
    # Assign max and min values
    vmin = 18.0
    vmax = 25.0
    # if flv == "MK10":
    #    vmin = 22.0
    #    vmax = 25.0
    # plot
    fig, ax = plt.subplots()
    h = ax.pcolormesh(TTs, PPs, np.log10(etapps), shading='auto', cmap='viridis_r',\
    vmin=vmin, vmax=vmax)
    ax.set_title("Peierls creep of %s, with eta = %.4e" % (flv, edotp0))
    ax.set_xlabel('Temperature (T)')
    ax.set_ylabel('Pressure (Pa)')
    fig.colorbar(h, ax=ax, label='log(Viscosity (Pa *s))')
    filename = "./Peierls_PT_%s_%s_edot%.4e.png" % (flv, implementation, edotp0)
    fig.savefig(filename)
    print("%s: export figure %s" % (Utilities.func_name(), filename))


def PlotPeierlsDFM():
    # Unit conversions
    mpa = 1e6  # MPa to Pa
    mum = 1e6  # microns/meter for grain size 
    R = 8.314       # J/mol*K
    km2m = 1e3 # km to meters
    sec2yrs = 60*60*24*365.25 # sec per year
    
    # Depth in meters
    zmin = 0
    zmax = 660
    
    depth = np.array([zmin, zmax])
    
    # Temperature range 
    mad = 0.5715  # deg/km, approximating adiabat to match T(CMB) in Aspect
    Tmin = 200
    Tmax = 2000
    dT = 50 # C
    
    Pmin = 0
    Pmax = 22 # GPa
    dP = 1 # GPa
    
    T1 = np.arange(Tmin,Tmax+dT,dT)
    P1 = np.arange(Pmin,Pmax+dP,dP)*1e9  # Pa
    
    T, P = np.meshgrid(T1,P1)
    
    # This is just here so I could look at what the pressure dependence gives for a change
    # in the peierls stress: goes from 5.9 to 8.6 GPa
    # From Kawazoe et al PEPI 2009
    G0 = 77.4*1e9 # GPa  
    Gp = 1.61 # GPa/GPa 
    sigp0 = 5.9e9
    sigp = sigp0*(1 + (Gp/G0)*P1)
    	
    # constant grain size and strain rate
    d = 10e3 	  # microns (= 10 mm = 1.0 cm)
    edot = 1e-15  # background mantle strain-Rate
    edot1 = 1e-15  # background mantle strain-Rate
    
    # Get viscosity
    water = 'wet'
    flver = 'mod'
    coh = 1000  
    
    # Diffusion Creep
    A, Emid1, Vmid1, p, r = get_diff_HK_params(water,flver,'mid','mid')    
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
    
    dE = -50e3      # dE = +/- 75, error on activation energy, from HK03
    dV = -4.5e-6  	# dV = +/-4.5e-6  # Ohuchi et al, 2012	
    
    Edf = Emid1 + dE
    Vdf = Vmid1 + dV	
    fh2o = convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etadf = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((Edf + P*Vdf)/(R*T)) # Pa s
    
    # Dislocation Creep
    A, Emid2, Vmid2, n, r = get_disl_HK_params(water,flver,'mid','mid')    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
    
    dE = -25e3    # dE = +/-40e3  error on activation energy from HK03
    dV = 3e-6  	 # dV = +/-3e-6 Ohuchi et al, 2012
    lab1 = 'df: ' + str(Edf/1e3) + '/' + str(Vdf*1e6)
    
    Eds = Emid2 + dE
    Vds = Vmid2 + dV
    lab2 = 'ds: ' + str(Eds/1e3) + '/' + str(Vds*1e6)	
    					
    fh2o = convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etads =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((Eds + P*Vds)/(n*R*T)) # Pa s
    	
    # Peierls Creep
    etap = peierls_approx_visc('MK10',P,T,edot)
    
    
    # Composite viscosity
    etacomp1 = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)
    etacomp = etadf*etads/(etadf+ etads)
    sig = 2*etacomp*edot
    p = np.where(sig>100e6)
    #etacomp[p] = etacomp1[p]
    etacomp = etadf*etads*etap/(etadf*etads + etadf*etap + etads*etap)

    # subfig1: diffusion creep viscosity	
    vmin1 = 18
    vmax1 = 25		
    fig = plt.figure(figsize=(8.5,11.0))
    ax1 = fig.add_subplot(3,2,1)	
    c1 = ax1.pcolor(T,P/1e9,np.log10(etadf),cmap='viridis')
    c1.cmap.set_over('k')
    c1.set_clim(vmin1,vmax1)
    ax1.set_ylabel('Pressure (GPa)')
    ax1.set_title('Diffusion',fontsize=10)	
    fig.colorbar(c1,ax=ax1)

    # subfig2: dislocation creep viscosity, the strain rate is T-dependent
    ax2 = fig.add_subplot(3,2,2)	
    c2 = ax2.pcolor(T,P/1e9,np.log10(etads),cmap='viridis')
    c2.cmap.set_over('k')
    c2.set_clim(vmin1,vmax1)
    ax2.set_title('Dislocation',fontsize=10)
    fig.colorbar(c2,ax=ax2)
 
    # subfig3: peierls creep, strain rate is T-dependent
    ax3 = fig.add_subplot(3,2,3)	
    c3 = ax3.pcolor(T,P/1e9,np.log10(etap), cmap='viridis')
    # ax3.contour(T,P/1e9,np.log10(etap), [21.0, 22.0], colors='white')
    ax3.contour(T,P/1e9,np.log10(etap), [23.0], colors='white')
    c3.cmap.set_over('k')
    c3.set_clim(vmin1,vmax1)
    ax3.set_xlabel('Temperature (C)')
    ax3.set_title('Peierls',fontsize=10)	
    fig.colorbar(c3,ax=ax3)

    # subfig 4: composite viscosity 
    ax4 = fig.add_subplot(3,2,4)	
    c4 = ax4.pcolor(T,P/1e9,np.log10(etacomp),cmap='viridis')
    c4.cmap.set_over('k')
    c4.set_clim(vmin1,vmax1)
    ax4.contour(T,P/1e9,np.log10(etacomp), [21.0, 22.0],colors='white')
    ax4.set_xlabel('Temperature (C)')
    ax4.set_title('Composite',fontsize=10)
    fig.colorbar(c4,ax=ax4)

    # subfig5: stress, with the strain rate dependent on T and the computed composite viscosity
    ax5 = fig.add_subplot(3,2,5)	
    c5 = ax5.pcolor(T,P/1e9,sig/1e9,cmap='viridis')
    c5.cmap.set_over('k')
    c5.set_clim(0,0.1)
    ax5.set_xlabel('Temperature (C)')
    ax5.set_title('Stress df-ds only',fontsize=10)
    fig.colorbar(c5,ax=ax5)
    
    skp = 1
    if skp == 0:
    	ax6 = fig.add_subplot(3,2,6)	
    	c6 = ax6.pcolor(T,P/1e9,np.log10(edot),cmap='viridis')
    	ax6.set_xlabel('Temperature (C)')
    	ax6.set_title('Strain-rate (T)',fontsize=10)
    	fig.colorbar(c6,ax=ax6)
    
    plt.tight_layout()
    pdffile = 'peierls_defmech_1e-15.pdf'
    fig.savefig(pdffile,bbox_inches='tight')
    print('%s: file generated %s' % (Utilities.func_name(), pdffile))
    assert(os.path.isfile(pdffile))



    
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
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    parser.add_argument('-im', '--implementation', type=str,
                        default='exact',
                        help='Implementation used for peierls creep')
    parser.add_argument('-f', '--flavor', type=str,
                        default='orig',
                        help='The flavor to use in computing of the HK03 rhoelogy')
    # parse options 
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)
    
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'plot_upper_mantle_viscosity':
        plot_upper_mantle_viscosity(ptplt=True, flavor=arg.flavor)
    elif _commend == 'plot_peierls_creep_PT': 
        edotp0 = 1e-15
        PlotPeierlsPTMap(arg.flavor, arg.implementation, edotp0)
    elif _commend == 'plot_peierls_dfm':
        PlotPeierlsDFM()
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)



# run script
if __name__ == '__main__':
    main()