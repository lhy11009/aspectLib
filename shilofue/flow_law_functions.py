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

# Peierls creep flow law 
# flv: flow law version
# MK10: for Mei and Kohlstedt, 2010	 (gam = 0.17)
# gam: fitting parameter = sig_ref/sigp
#
def peierls_approx_visc(flv,gam,P,T,edot):

	mpa = 1e6  # MPa to Pa

	if flv == "MK10":
		# Mei et al., JGR 2010
		q = 1.0
		p = 0.5
		n = 2.0
		sigp0 = 5.9e9					# Pa (+/- 0.2e9 Pa)
		A = 1.4e-7/np.power(mpa,n) 	# s^-1 Pa^-2
		E = 320e3  					# J/mol (+/-50e3 J/mol)
	
	# Pressure dependence of Peierls stress
	# From Kawazoe et al. PEPI 2009 and parameters from Liu et al., GRL 2005 
	G0 = 77.4*1e9 # GPa  
	Gp = 1.61 # GPa/GPa 
	sigp = sigp0*(1 + (Gp/G0)*P)
	
	s = (E/(R*T))*p*q*((1-gam**p)**(q-1))*(gam**p)
	x = 1/(s+n)
	visc = (0.5*gam*sigp*edot**(x-1))/( ((A*(gam*sigp)**n)**x)*np.exp( (E*(1-gam**p)**q)/(R*T*(s+n)) ) )
	
	return visc

    
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
    parser.add_argument('-f', '--flavor', type=str,
                        default='orig',
                        help='The flavor to use in computing of the HK03 rhoelogy')
    # parse options 
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)
    
    if _commend == 'plot_upper_mantle_viscosity':
        plot_upper_mantle_viscosity(ptplt=True, flavor=arg.flavor)



# run script
if __name__ == '__main__':
    main()