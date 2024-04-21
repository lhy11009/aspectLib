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