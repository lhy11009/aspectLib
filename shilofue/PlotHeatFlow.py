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

# todo
def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - plot the heat flow of the upper bd: \n\
\n\
        python -m \
        ")


class BOUNDARYOUTPUTs():
    '''
    A class that handles boundary data output, including reading data files, processing boundaries based on geometry,
    and exporting heat flow data for specified boundaries.

    Attributes:
    ndim (int): Dimension of the data, typically 2 for 2D.
    xs (np.ndarray or None): Array to store x-coordinates of boundary points.
    ys (np.ndarray or None): Array to store y-coordinates of boundary points.
    zs (np.ndarray or None): Array to store z-coordinates of boundary points (set to zero for 2D).
    Rs (np.ndarray or None): Array to store radial distances for boundary points.
    Thetas (np.ndarray or None): Array to store theta angles for spherical coordinates.
    Phis (np.ndarray or None): Array to store phi angles (longitude) for spherical coordinates.
    hfs (np.ndarray or None): Array to store heat flux values at the boundary points.
    geometry (str): Type of geometry for the domain, either "box" or "chunk" (default is "box").
    bd_indicators (np.ndarray or None): Array that holds boundary indicators (0 for left, 1 for right, 2 for bottom, 3 for top).
    '''
    
    def __init__(self, ndim: int, **kwargs):
        '''
        Initializes the object with default values for coordinates, heat flux, and geometry.
        
        Parameters:
        ndim (int): The dimension of the data (2D by default).
        kwargs (dict): Additional keyword arguments for geometry options (default is "box").
        '''
        self.ndim = ndim
        self.xs = None
        self.ys = None
        self.zs = None
        self.Rs = None
        self.Thetas = None
        self.Phis = None
        self.hfs = None
        self.geometry = kwargs.get("geometry", "box")
        self.bd_indicators = None

    def ReadFile(self, filein: str):
        '''
        Reads the boundary data from a file. The file is expected to contain columns of x, y, and heat flux (hfs) data.
        For 2D data, z is initialized as 0. Ensures the file exists and loads the data into the class attributes.
        
        Parameters:
        filein (str): The input file path from which to read boundary data.
        '''
        Utilities.my_assert(os.path.isfile(filein), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), filein))
        data = np.genfromtxt(filein)
        if self.ndim == 2:
            self.xs = data[:, 0]
            self.ys = data[:, 1]
            self.zs = np.zeros(self.xs.shape)
            self.hfs = data[:, 2]
        else:
            raise NotImplementedError

    def ProcessBds(self, **kwargs):
        '''
        Processes the boundary data depending on the geometry type. Ensures that heat flux data (hfs) exists. 
        If the geometry is "chunk", it processes the boundaries for this geometry; otherwise, it raises an error.
        
        Parameters:
        kwargs (dict): Additional parameters for processing boundaries (e.g., inner_radius, outer_radius).
        '''
        Utilities.my_assert(self.hfs is not None, ValueError, "%s: no data present." % (Utilities.func_name()))
        if self.geometry == "chunk":
            self.ProcessBdsChunk(**kwargs)
        else:
            raise NotImplementedError
    
    def ProcessBdsChunk(self, **kwargs):
        '''
        Processes boundaries for "chunk" geometry by assigning boundary indicators and converting Cartesian coordinates 
        to spherical coordinates. It calculates radial distances and assigns boundary indicators for specific boundaries 
        based on radius (Ri, Ro) and longitude (phi0, phi1).
        
        Parameters:
        kwargs (dict): Contains inner_radius (float), outer_radius (float), minimum_longitude (float),
                       and maximum_longitude (float) for boundary processing.
        '''
        Ri = kwargs.get("inner_radius", 3.481e6)
        Ro = kwargs.get("outer_radius", 6.371e6)
        phi0 = kwargs.get("minimum_longitude", 0.0)
        phi1 = kwargs.get("maximum_longitude", 1.4002e+02)
        
        self.Rs = (self.xs**2.0 + self.ys**2.0)**0.5
        self.Rs, self.Thetas, self.Phis = Utilities.cart2sph(self.xs, self.ys, self.zs)

        self.bd_indicators = np.full(self.Rs.shape, np.nan)

        # Assign boundary indicators based on longitude and radius thresholds.
        # 0 - left, 1 - right, 2 - bottom, 3 - top boundaries.
        if self.ndim == 2:
            mask0 = (np.abs(self.Phis - phi0) < 1e-3)
            self.bd_indicators[mask0] = 0
            mask1 = (np.abs(self.Phis - phi1) < 1e-3)
            self.bd_indicators[mask1] = 1
            mask2 = (np.abs(self.Rs - Ri) < 1e3)
            self.bd_indicators[mask2] = 2
            mask3 = (np.abs(self.Rs - Ro) < 1e3)
            self.bd_indicators[mask3] = 3
        else:
            raise NotImplementedError

    def ExportBdHeatFlow(self, indicator: int):
        '''
        Exports heat flow data for a specified boundary based on the boundary indicator. It returns the x, y, z coordinates
        and heat flux values for the boundary specified by the indicator.
        
        Parameters:
        indicator (int): The boundary indicator (0 - left, 1 - right, 2 - bottom, 3 - top).
        
        Returns:
        Tuple: Contains arrays of x, y, z coordinates and corresponding heat flux values.
        '''
        mask = (self.bd_indicators == indicator)
        return self.xs[mask], self.ys[mask], self.zs[mask], self.hfs[mask]


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
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()