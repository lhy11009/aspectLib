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
import math
# import json, re
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

def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
    pass

class HELE_SHAW():

    def __init__(self, Vt, Vm, Vr, L, lbd, mu, n):
        '''
        initiation
        Attributes:
            Vt - upper plate velocity 
            Vm - lower mantle velocity
            Vr - trench velocity
            L - length of the slab 
            lbd - thickness of the upper mantle
            mu - viscosity
            n - number of points along half the slab
        '''
        self.Vt = Vt
        self.Vm = Vm
        self.Vr = Vr
        self.a = L / 2.0
        self.lbd = lbd
        self.mu = mu
        self.n = n
        self.As = np.zeros(self.n)
        # self.As = 0.001 * np.ones(self.n) * mu / lbd**2.0 * (Vt + Vm) # for test
        pass

    def AssignA(self, As):
        '''
        Assign an initial guess of As
        Inputs:
            As - an array of summation coefficiences
        '''
        self.As = As

    def Fxiyi(self, x, y, xi, yi):
        '''
        coefficients of the solution for Vx, (partial P, partial x)
        Inputs:
            x - coordinate x
            y - coordinate y
            xi - coordinate x of i th summation point
            yi - coordinate y of i th summation point
        '''
        temp = ((x - xi)**2.0 - (y - yi)**2.0) / ((x - xi)**2.0 + (y - yi)**2.0)**2.0
        return temp
    
    def F1xiyi(self, x, y, xi, yi):
        '''
        coefficients of the solution for Vy, (partial P, partial y)
        Inputs:
            x - coordinate x
            y - coordinate y
            xi - coordinate x of i th summation point
            yi - coordinate y of i th summation point
        '''
        temp = (2.0*(x - xi)*(y - yi)) / ((x - xi)**2.0 + (y - yi)**2.0)**2.0
        return temp
    
    def Gxiyi(self, x, y, xi, yi):
        '''
        coefficients of the solution for P
        Inputs:
            x - coordinate x
            y - coordinate y
            xi - coordinate x of i th summation point
            yi - coordinate y of i th summation point
        '''
        temp = (x - xi) / ((x - xi)**2.0 + (y - yi)**2.0)
        return temp
    
    def Vx(self, x, y, z):
        '''
        Solution for Vx
        '''
        temp = 0.0
        for i in range(self.n):
            Ai = self.As[i]
            xi = 0.0
            yi = i / (self.n-1.0) * self.a
            yi_n = -i / (self.n-1.0) * self.a
            temp += (self.Fxiyi(x, y, xi, yi) if (i == 0) 
                    else self.Fxiyi(x, y, xi, yi_n) + self.Fxiyi(x, y, xi, yi))\
                    * Ai / (2.0 * self.mu)
        Vx = (self.Vt * (z - self.lbd) + self.Vm*z) / self.lbd + temp * z * (self.lbd - z)
        return Vx
    
    def Vy(self, x, y, z):
        '''
        Solution for Vx
        '''
        temp = 0.0
        for i in range(self.n):
            Ai = self.As[i]
            xi = 0.0
            yi = i / (self.n-1.0) * self.a
            yi_n = -i / (self.n-1.0) * self.a
            temp += (self.F1xiyi(x, y, xi, yi) if (i == 0) 
                    else self.F1xiyi(x, y, xi, yi_n) + self.F1xiyi(x, y, xi, yi))\
                    * Ai / (2.0 * self.mu)
        Vy = temp * z * (self.lbd - z)
        return Vy


    def Vxdz(self, x, y):
        '''
        integration of Vx over z
        Inputs:
            x - coordinate x
            y - coordinate y
        '''
        sum = self.lbd * (self.Vt + self.Vm) / 2.0
        for i in range(self.n):
            Ai = self.As[i]
            xi = 0.0
            yi = i / (self.n-1.0) * self.a
            yi_n = -i / (self.n-1.0) * self.a
            # sum over the positive value yi and negative value yi_n
            if x == np.inf:
                sum += (1.0 if (i == 0) else 2.0)\
                        * self.lbd**3.0 * Ai / (12.0 * self.mu)
            else:
                sum += (self.Fxiyi(x, y, xi, yi) if (i == 0) 
                        else self.Fxiyi(x, y, xi, yi_n) + self.Fxiyi(x, y, xi, yi))\
                      * self.lbd**3.0 * Ai / (12.0 * self.mu)
        return sum
    
    def Pressure(self, x, y):
        '''
        integration of Vx over z
        Inputs:
            x - coordinate x
            y - coordinate y
        '''
        sum = self.lbd * (self.Vt + self.Vm) / 2.0
        for i in range(self.n):
            Ai = self.As[i]
            xi = 0.0
            yi = i / (self.n-1.0) * self.a
            yi_n = -i / (self.n-1.0) * self.a
            # sum over the positive value yi and negative value yi_n
            sum += (self.Gxiyi(x, y, xi, yi) if (i == 0) 
                    else self.Gxiyi(x, y, xi, yi_n) + self.Gxiyi(x, y, xi, yi))\
                      * Ai
        return sum
    
    def Vxdz_array(self, x):
        '''
        derive an array of Vxdz over different value of yi
        Inputs:
            x - coordinate x
        '''
        Vxdzs = np.zeros(self.n)
        for i in range(self.n):
            y = i * self.a / (self.n - 1.0)
            Vxdzs[i] = self.Vxdz(x, y)
        return Vxdzs

    def Vxdz_array_nd(self, x):
        '''
        derive the nondimentional array of Vxdz
        Inputs:
            x - coordinate x
        '''
        return self.Vxdz_array(x) / (self.Vr * self.lbd)
    
    def Vxdz_array_nd_full(self, x, **kwargs):
        '''
        derive the nondimentional array of Vxdz, include the value at infinite x
        Inputs:
            x - coordinate x
            kwargs:
                weight_inf - weight of the infinite point
        '''
        weight_inf = kwargs.get("weight_inf", 0.0)

        temp = np.zeros(self.n+1)
        temp[0:self.n] = self.Vxdz_array_nd(x)
        temp[self.n] = weight_inf * self.Vxdz(np.inf, 0.0) / (self.Vr * self.lbd)
        return temp
    
    def partial_Vxdz_i_partial_Aj_nd(self, x, i, j):
        '''
        partial derivative of Vxdz at the ith point to Aj
        Inputs:
            x - coordinate x
            i - index i
            j - index j
        '''
        temp = (j - i) * (j - i) * self.a * self.a / (self.n-1) / (self.n-1)
        return (x * x - temp) / (x * x + temp)**2.0 * self.lbd * self.lbd / (12 * self.mu * self.Vr)
    
    def partial_Vxdz_i_partial_Aj_nd_xinf(self):
        '''
        partial derivative of Vxdz at infinite x
        '''
        return self.lbd * self.lbd / (12 * self.mu * self.Vr)

    def Jacobi_slab_ij(self, x, i, j):
        '''
        Jacobian of partial Vxdzs_i partial Aj, only include the points at the slab surface
        '''
        J = (self.partial_Vxdz_i_partial_Aj_nd(x, i, j) if (j == 0)\
            else self.partial_Vxdz_i_partial_Aj_nd(x, i, j) + self.partial_Vxdz_i_partial_Aj_nd(x, -i, j))
        return J
    
    def Jacobi_slab(self, x):
        '''
        Full Jacobian matrix of partial Vxdzs_i partial Aj, only include the points at the slab surface
        Inputs:
            x - coordinate x
        '''
        J = np.zeros([self.n, self.n])
        # sum over the positive indexes i and negative value -i
        for i in range(self.n):
            for j in range(self.n):
                J[i, j] = self.Jacobi_slab_ij(x, i, j)
        return J
    
    def Jacobi_full(self, x, **kwargs):
        '''
        Jocobian of partial Vxdzs_i partial Aj, include also the infinite distance along x
        Inputs:
            x - coordinate x
            kwargs:
                weight_inf - weight of the infinite point
        '''
        weight_inf = kwargs.get("weight_inf", 0.0)
        
        J = np.zeros([self.n+1, self.n])
        # sum over the positive indexes i and negative value -i
        for i in range(self.n):
            for j in range(self.n):
                J[i, j] = self.Jacobi_slab_ij(x, i, j)
        # additional component at x = inf
        for j in range(self.n): 
            J[self.n, j] = weight_inf * (1.0 if (j==0) else 2.0) * self.partial_Vxdz_i_partial_Aj_nd_xinf()
        return J
    
    def PlotVxdz(self, x, ym, **kwargs):
        '''
        Plot the solution
        Inputs:
            x - coordinate x
            ym - maximum value of y
            kwargs:
                ax - plotting axis
        '''
        ax = kwargs.get("ax", None)
        N = 1000
        ys = np.linspace(0, ym, N)
        Vxdzs = np.zeros(N)
        for i in range(N):
            y = ys[i]
            Vxdzs[i] = self.Vxdz(x, y)
            Vx_dzs_nd = Vxdzs / self.lbd / self.Vr
        if ax is not None:
            ax.plot(ys / self.a, Vxdzs / self.lbd / self.Vr)
        ax.set_xlabel("Distance from center of trench (y / a)")
        ax.set_ylabel("Velocity at x=0 (Vx/Vr)")

    def PlotCenterPressureX(self, xm, **kwargs):
        '''
        Plot the solution
        Inputs:
            x - coordinate x
            ym - maximum value of y
            kwargs:
                ax - plotting axis
        '''
        ax = kwargs.get("ax", None)
        x_nd = kwargs.get('x_nd', 10e3)
        y_nd = 0.0

        # pressure value used for non-dimentionlization
        P0 = self.Pressure(x_nd, y_nd)

        N = 1000
        xs = np.linspace(0, xm, N)
        Ps = np.zeros(N)
        for i in range(N):
            x = xs[i]
            y = 0.0
            Ps[i] = self.Pressure(x, y)
        if ax is not None:
            ax.plot(xs / self.a, Ps / P0)
        ax.set_xlabel("Distance from slab (x / a)")
        ax.set_ylabel("Normalized Pressure at y=0")
    
    def PlotVelocityOnSurface(self, z0, **kwargs):
        '''
        plot the velocity on the surface at depth z0
        '''
        yr = 365.0 * 24.0 * 3600.0
        cm_per_yr = 1e-2 / yr
        vector_scale = 1.0 / cm_per_yr

        ax = kwargs.get("ax", None)
        xmin = kwargs.get("xmin", 10e3)
        xmax = kwargs.get("xmax", 2000e3)
        nx = kwargs.get("nx", 20)
        ymin = kwargs.get("ym", 0e3)
        ymax = kwargs.get("ym", 3.0 * self.a)
        ny = kwargs.get("ny", 15)

        xs, ys = symmetric_grid(xmin, xmax, nx, ymin, ymax, ny)

        xxs, yys = np.meshgrid(xs, ys)
        vvxs = np.zeros(xxs.shape)
        vvys = np.zeros(yys.shape)

        for i in range(ys.size):
            for j in range(xs.size):
                vvxs[i, j] = self.Vx(xs[j], ys[i], z0)
                vvys[i, j] = self.Vy(xs[j], ys[i], z0)

        Q = ax.quiver(xxs/1e3, yys/1e3, vvxs * vector_scale, vvys * vector_scale, angles='xy')
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('X [km]')
        ax.set_ylabel('Y [km]')
        qk = ax.quiverkey(Q, 0.7, 0.92, 1.0, r'$1 \frac{cm}{yr}$', labelpos='E',
                   coordinates='figure')
        return xxs, yys, vvxs * vector_scale, vvys * vector_scale
    
    def PlotPressureOnSurface(self, **kwargs):
        '''
        plot the velocity on the surface at depth z0
        '''
        Mpa = 1e6

        ax = kwargs.get("ax", None)
        xmin = kwargs.get("xmin", 10e3)
        xmax = kwargs.get("xmax", 2000e3)
        nx = kwargs.get("nx", 20)
        ymin = kwargs.get("ym", 0e3)
        ymax = kwargs.get("ym", 3.0 * self.a)
        ny = kwargs.get("ny", 15)
        vmin = kwargs.get("vmin", None)
        vmax = kwargs.get("vmax", None)

        xs, ys = symmetric_grid(xmin, xmax, nx, ymin, ymax, ny)

        xxs, yys = np.meshgrid(xs, ys)
        PPs = np.zeros(xxs.shape)

        for i in range(ys.size):
            for j in range(xs.size):
                PPs[i, j] = self.Pressure(xs[j], ys[i])

        h = ax.pcolormesh(xxs/1e3, yys/1e3, PPs/Mpa, vmin=vmin, vmax=vmax)
        return h

def symmetric_grid(xmin, xmax, nx, ymin, ymax, ny):
    '''
    Generate symmeric grid on a plane of a contant z
    '''
    xs1 = np.linspace(xmin, xmax, nx)
    xs2 = np.linspace(-xmax, xmin, nx)
    xs = np.concatenate((xs1, xs2))
    ys = np.linspace(-ymax, ymax, 3*ny+1)
    return xs, ys

class SLAB():
    '''
    class for slab profile
    Attributes:
        Vs: velocity along s
        Vn: velocity along n
        Taus: shear stress
        SigNs: normal stress
        pSigN_pSs: normal stress derivatives
        P0: viscous pressure from large scale flow around the slab
    '''
    def __init__(self, lbd, mu):
        self.n = 100
        self.lbd = lbd
        # initiate the profile of x and y
        self.xs = np.zeros(self.n)
        self.zs = np.linspace(0, self.lbd, self.n)
        self.r0s = np.zeros(self.n)
        self.theta0s = np.zeros(self.n)
        self.Vs = 0.0
        self.Vns = np.zeros(self.n)
        self.pSigN_pSs = np.zeros(self.n)
        self.Taus = np.zeros(self.n)
        self.mu = mu
        self.SigNs = np.zeros(self.n)
        self.P0 = 0.0

    def LocalProfile(self, i):
        '''
        '''
        assert(i < self.n)
        # value of the interpolation points
        x0 = 0.0
        y0 = 0.0
        z0 = 0.0
        z1 = 0.0
        if i == 0:
            x0 = self.xs[0]
            x1 = self.xs[1]
            z0 = self.zs[0]
            z1 = self.zs[1]
        elif i == self.n-1:
            x0 = self.xs[self.n-2]
            x1 = self.xs[self.n-1]
            z0 = self.zs[self.n-2]
            z1 = self.zs[self.n-1]
        else:
            x0 = self.xs[i-1]
            x1 = self.xs[i+1]
            z0 = self.zs[i-1]
            z1 = self.zs[i+1]
        theta0 = math.atan2(x1 -x0, z1-z0)
        r0 = (self.lbd - self.zs[i]) / math.asin(theta0)
        return r0, theta0

    def UpdateProfile(self):
        '''
        Update the slab profile
        '''
        for i in range(self.n):
            r0, theta0 = self.LocalProfile(i)
            self.r0s[i] = r0
            self.theta0s[i] = theta0
    
    def PlotProfile(self, **kwargs):
        '''
        plot the slab profile
        '''
        ax = kwargs.get('ax', None)
        ax.plot(self.xs/1e3, self.zs/1e3)
        # plot properties
        ax.set_xlabel("X (km)")
        ax.set_xlim([0.0, self.lbd*2.0/1e3])
        ax.set_ylabel("Z (km)")
        ax.set_ylim([0.0, self.lbd/1e3])
        ax.invert_xaxis()
        ax.invert_yaxis()
        ax.grid()

    def InitiateProfile(self, S, theta):
        '''
        assign an initial profile
        '''
        for i in range(self.n):
            self.xs[i] = S * i / (self.n - 1.0) * math.cos(theta)
            self.zs[i] = S * i/ (self.n - 1.0) * math.sin(theta)
        self.UpdateProfile()
    
    def UpdateVs(self, Vs):
        '''
        update Vs on the slab
        '''
        self.Vs = Vs

    def calculate_SigN_Tau(self):
        '''
        eq 5, 6
        '''
        # sum over Vn
        sum_Vns = np.zeros(self.n)
        # derive normal and shear stress
        for i in range(self.n):
            r0 = self.r0s[i]
            theta0 = self.theta0s[i]
            Vn = self.Vns[i]
            Vs = self.Vs
            a = -12.0 * self.mu / (r0**3.0 * theta0**3.0) * sum_Vns[i]
            b = 6 * self.mu * (Vn + Vs) / r0**2.0 / theta0**2.0
            self.pSigN_pSs[i] = a + b  
            self.Taus[i] = -r0 * theta0 * self.pSigN_pSs[i] + self.mu * (Vn - Vs) / r0 / theta0
        self.SigNs[0] = self.P0
        for i in range(1, self.n):
            dz = self.zs[i] - self.zs[i-1]
            self.SigNs[i] = self.SigNs[i-1] + (self.pSigN_pSs[i-1] + self.pSigN_pSs[i]) * dz / 2.0
    
    def calculate_P0(self):
        '''
        eq 4
        '''
        pass

    def PlotStresses(self, **kwargs):
        '''
        Plot the stress profile
        '''
        ax = kwargs.get('ax', None)
        ax.plot(self.SigNs / 1e6, self.zs/1e3, label="normal")
        ax.plot(self.Taus / 1e6, self.zs/1e3, label="shear")
        # plot properties
        ax.set_xlabel("Stress (Mpa)")
        ax.set_ylabel("Z (km)")
        ax.set_ylim([0.0, self.lbd/1e3])
        ax.invert_yaxis()
        ax.legend()
        ax.grid()


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