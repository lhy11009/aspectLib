# -*- coding: utf-8 -*-
#!/opt/anaconda3/bin/python3
r"""Test for peierls rheology

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest 

descriptions:
    every function is a separate test, combined usage with pytest module
""" 

# makes  a deformation mechanism map in P-T space for constant grain size and strainrate

import os
import shutil
import math
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import shilofue.FlowLaws as flf

test_dir = ".test"
# source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_run_time')

def test_get_peierls():
    '''
    test the interface for exporting peierls creep
    '''
    Peierls = flf.GetPeierlsApproxVist('MK10')
    assert(abs(Peierls['q'] - 1.0) < 1e-6)
    assert(abs((Peierls['sigp0'] - 5900000000.0)/5900000000.0) < 1e-6)
    assert(abs((Peierls['A'] - 1.4e-19)/1.4e-19) < 1e-6)
    # print(Peierls)  # print what's there


def test_peierls_visc():
    '''
    test peierls_visc_from_stress function
    assert:
        values of strain rate
    '''
    # test 1
    P = 0.0
    T = 1050.0
    sigma = 3e9
    etap, edotp = flf.peierls_visc_from_stress('MK10',P,T, sigma)
    assert(abs((edotp - 3.408743e-05)/3.408743e-05) < 1e-6)
    # test 2
    P = 0.0
    T = 1400.0
    sigma = 2.35e9
    etap, edotp = flf.peierls_visc_from_stress('MK10',P,T, sigma)
    assert(abs((edotp - 3.0468359e-05)/3.0468359e-05) < 1e-6)


def test_peierls_visc_from_edot():
    '''
    test peierls_visc_from_edot function
    assert:
        values of viscosity and n (number of iteration)
    '''
    # test 1
    P = 0.0
    T = 1050.0
    edot = 3.408743e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot('MK10',P,T, edot)
    assert(abs((etap - 4.480804e+13)/4.480804e+13) < 1e-5)
    assert(n == 7)
    # test 2
    P = 0.0
    T = 1400.0
    edot = 3.0468359e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot('MK10',P,T, edot)
    assert(abs((etap - 4.069405e+13)/4.069405e+13) < 1e-5)
    assert(n == 11)


def test_peierls_visc_from_edot_newton_nolog():
    '''
    test peierls_visc_from_edot_newton_nolog function
    assert:
        values of viscosity and n (number of iteration)
    '''
    # test 0
    P = 1.7882e9
    T = 1330.38
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('MK10',P,T, edot, 0.1)
    assert(abs((etap - 2.167156389993528e+22)/2.167156389993528e+22) < 1e-5)
    assert(n == 5)
    # test 1
    P = 0.0
    T = 1050.0
    edot = 3.408743e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('MK10',P,T, edot, 0.1)
    assert(abs((etap - 4.26878809748e+13)/4.26878809748e+13) < 1e-5)
    assert(n == 6)
    # test 2
    P = 0.0
    T = 1400.0
    edot = 3.0468359e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('MK10',P,T, edot, 0.1)
    assert(abs((etap - 3.737963093140574e+13)/3.737963093140574e+13) < 1e-5)
    assert(n == 3)


def test_peierls_newton_log_nolog_compare():
    '''
    use the peierls_visc_from_edot_newton function 
    and the peierls_visc_from_edot_newton_nolog function
    to compare the use of log value and value as gradients
    '''
    # test 1.1: nolog, medium T
    P = 1.7882e9
    T = 1330.38
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('MK10',P,T, edot, 1e-8)
    assert(abs((etap - 2.273489862826893e+22)/2.273489862826893e+22) < 1e-5)
    assert(n == 8)
    # test 1.2: log, medium T
    P = 1.7882e9
    T = 1330.38
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton('MK10',P,T, edot, 1e-8)
    assert(abs((etap - 2.2734917295040908e+22)/2.2734917295040908e+22) < 1e-5)
    assert(n == 4)
    pass
    # test 2.1: nolog, high T
    P = 1.7882e9
    T = 3500.0
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('MK10',P,T, edot, 1e-8)
    assert(abs((etap - 1.0219201878630164e+19)/1.0219201878630164e+19) < 1e-5)
    assert(n == 11)
    # test 2.2: log, high T
    P = 1.7882e9
    T = 3500.0
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton('MK10',P,T, edot, 1e-8)
    assert(abs((etap - 1.0219201878630164e+19)/1.0219201878630164e+19) < 1e-5)
    assert(n == 3)
    # test 3.1: nolog, high T, Idrissi
    P = 1.7882e9
    T = 3500.0
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton_nolog('Idrissi16',P,T, edot, 1e-8)
    assert(math.isnan(etap))
    assert(n == 1)
    # test 3.2: log, high T, Idrissi
    # doesn't converge, returns to sigma = 0.0 and edotp != 0.0
    P = 1.7882e9
    T = 3500.0
    edot = 1e-15
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton('Idrissi16',P,T, edot, 1e-8, debug=True)
    assert(abs(etap) < 1e-5)


def test_peierls_visc_from_edot_newton():
    '''
    test peierls_visc_from_edot_newton function
    assert:
        values of viscosity and n (number of iteration)
    '''
    # test 1
    P = 0.0
    T = 1050.0
    edot = 3.408743e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton('MK10',P,T, edot, 0.1)
    assert(abs((etap - 4.396405536e+13)/4.396405536e+13) < 1e-5)
    assert(n == 2)
    # test 2
    P = 0.0
    T = 1400.0
    edot = 3.0468359e-05
    etap, sigma, diff, n = flf.peierls_visc_from_edot_newton('MK10',P,T, edot, 0.1)
    assert(abs((etap - 3.67999922e+13)/3.67999922e+13) < 1e-5)
    assert(n == 1)


def test_plot_peierls():
    # Unit conversions
    mpa = 1e6  # MPa to Pa


    peierls_visc_from_edot_newton


def test_plot_peierls():
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
    #edot = 1e-15  # background mantle strain-Rate
    edot = 10**((3/1000)*T - 19)  # haoyuan: are they dependent on temperature?
    edot1 = 10**((3/1000)*T1 - 19)
    
    # Get viscosity
    water = 'wet'
    flver = 'mod'
    coh = 1000  
    
    # Diffusion Creep
    A, Emid1, Vmid1, p, r = flf.get_diff_HK_params(water,flver,'mid','mid')    
    dm = d/1e6  # convert from microns to meters
    Am = A/1e6  # convert from MPa^-1 to Pa^-1
    
    dE = -50e3      # dE = +/- 75, error on activation energy, from HK03
    dV = -4.5e-6  	# dV = +/-4.5e-6  # Ohuchi et al, 2012	
    
    Edf = Emid1 + dE
    Vdf = Vmid1 + dV
    print("Diffusion E: ", Edf)	
    print("Diffusion V: ", Vdf)	  # debug
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etadf = 0.5*(dm)**p/(Am*fh2o**r)*np.exp((Edf + P*Vdf)/(R*T)) # Pa s
    
    # Dislocation Creep
    A, Emid2, Vmid2, n, r = flf.get_disl_HK_params(water,flver,'mid','mid')    
    Am = A/(1e6)**n  # convert from MPa^-1 to Pa^-1
    
    dE = -25e3    # dE = +/-40e3  error on activation energy from HK03
    dV = 3e-6  	 # dV = +/-3e-6 Ohuchi et al, 2012
    lab1 = 'df: ' + str(Edf/1e3) + '/' + str(Vdf*1e6)
    
    Eds = Emid2 + dE
    Vds = Vmid2 + dV
    print("Dislocation E: ", Eds)	
    print("Dislocation V: ", Vds)	  # debug
    lab2 = 'ds: ' + str(Eds/1e3) + '/' + str(Vds*1e6)	
    					
    fh2o = flf.convert_OH_fugacity(T,P,coh)	# Convert coh to water fugacity
    etads =  0.5*edot**(1/n - 1)/(Am*fh2o**r)**(1/n)*np.exp((Eds + P*Vds)/(n*R*T)) # Pa s
    	
    # Peierls Creep
    etap = flf.peierls_approx_visc('MK10', P,T,edot)
    
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
    c3 = ax3.pcolor(T,P/1e9,np.log10(etap),cmap='viridis')
    ax3.contour(T,P/1e9,np.log10(etap),colors='white')
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
    ax4.contour(T,P/1e9,np.log10(etacomp),colors='white')
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
    pdffile = os.path.join(test_dir, 'peierls_defmech.pdf')
    if os.path.isfile(pdffile):
      os.remove(pdffile)  # remove old file
    fig.savefig(pdffile,bbox_inches='tight')
    assert(os.path.isfile(pdffile))
