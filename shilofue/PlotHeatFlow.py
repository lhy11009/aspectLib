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
import pandas as pd
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
# from matplotlib import cm
from scipy.interpolate import interp1d, UnivariateSpline
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


def HeatFlowRetriveProfile(local_dir: str, _time: float, Visit_Options: object, **kwargs: dict) -> tuple:
    '''
    Retrieves the heat flow profile based on the given parameters.

    Parameters:
        local_dir (str): Directory where output files are located.
        _time (float): The time at which to retrieve the heat flow profile.
        Visit_Options (object): An object containing options for retrieving data.
        kwargs (dict): Additional keyword arguments; accepts 'phi_diff' (float, default: 5.0).

    Returns:
        tuple: A tuple containing two masked arrays: heat fluxes and corresponding phi values.
    '''
    phi_diff = kwargs.get("phi_diff", 5.0)  # in degree

    _time1, timestep, _ = Visit_Options.get_timestep_by_time(_time)
    filein = os.path.join(local_dir, "output", "heat_flux.%05d" % timestep)
    Utilities.my_assert(os.path.isfile(filein), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), filein))
    slab_morph_path = os.path.join(local_dir, "vtk_outputs", "slab_morph_t1.00e+05.txt")
    Utilities.my_assert(os.path.isfile(slab_morph_path), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), slab_morph_path))
    print("Plot file %s, time %.2f, timestep %d" % (filein, _time1, timestep))

    # Retrieve trench location from slab morphology data, 
    # ensure the slab morphology file exists, load the data, 
    # and interpolate the trench position based on the specified time.
    data = np.loadtxt(slab_morph_path)
    steps = data[:, 1]
    times = data[:, 2]
    trenches = data[:, 3]
    slab_depths = data[:, 4]
    sfunc = interp1d(times, trenches, assume_sorted=True)
    trench = sfunc(_time1)

    # Interpret Visit_Options, retrieve geometric parameters (inner and outer radii),
    # and initialize BOUNDARYOUTPUTs to handle boundary data.
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    Ro = Visit_Options.options['OUTER_RADIUS']  # Outer radius
    Ri = Visit_Options.options['INNER_RADIUS']  # Inner radius
    BdOutputs = BOUNDARYOUTPUTs(2, geometry=geometry)

    # Read the boundary output data from the specified file, process it 
    # to compute heat flux between inner and outer radii, and export 
    # the heat flow data (xs, ys, zs: cartesian coordinates; hfs: heat flux).
    BdOutputs.ReadFile(filein)
    BdOutputs.ProcessBds(inner_radius=Ri, outer_radius=Ro)
    xs, ys, zs, hfs = BdOutputs.ExportBdHeatFlow(3)
    Rs, Thetas, Phis = Utilities.cart2sph(xs, ys, zs)

    # Let's get the masked data adjacent to the trench
    phi0 = trench - phi_diff / 180.0 * np.pi  # Trench lower boundary
    phi1 = trench + phi_diff / 180.0 * np.pi  # Trench upper boundary
    mask = ((Phis > phi0) & (Phis < phi1))
    hfs_masked = hfs[mask]
    Phis_masked = Phis[mask]

    return hfs_masked, Phis_masked


def HeatFlowRetriveZeroCrossings(Phis_masked: np.ndarray, hfs_spline: object, output_path: str) -> np.ndarray:
    '''
    Identifies the zero crossings of the heat flow spline and saves the results to a CSV file.

    Parameters:
        Phis_masked (np.ndarray): Masked phi values corresponding to the heat flow.
        hfs_spline (object): A spline object representing the heat flow values.
        output_path (str): The file path where the results will be saved.

    Returns:
        np.ndarray: An array of indices where zero crossings occur.
    '''
    # Create a dense grid for x to evaluate the spline
    Phis_dense = np.linspace(Phis_masked.min(), Phis_masked.max(), 1000)  # 1000 points between min and max of longitude
    hfs_dense_value = hfs_spline(Phis_dense)
    dhf_dphi = hfs_spline.derivative(n=1)
    dhf_dphi2 = hfs_spline.derivative(n=2)
    dhf_dphi_dense_value = dhf_dphi(Phis_dense)
    dhf_dphi2_dense_value = dhf_dphi2(Phis_dense)

    # Find indices where the second derivative is approximately zero
    zero_crossings = np.where(np.isclose(dhf_dphi_dense_value, 0, atol=1e-2))[0]

    # Create a DataFrame from the arrays
    # Then save the DataFrame to a CSV file
    df = pd.DataFrame({'zero crossings': zero_crossings, 'longitude': Phis_dense[zero_crossings], 
                       'heat flow': hfs_dense_value[zero_crossings],
                       'first derivative': dhf_dphi_dense_value[zero_crossings], 
                       'second derivative': dhf_dphi2_dense_value[zero_crossings]})
    if not os.path.isdir(os.path.basename(output_path)):
        os.mkdir(os.path.basename(output_path))
    df.to_csv(output_path, index=False)
    print(f"Output heat flux csv file: {output_path}")

    return [zero_crossings, Phis_dense[zero_crossings], hfs_dense_value[zero_crossings],\
        dhf_dphi_dense_value[zero_crossings], dhf_dphi2_dense_value[zero_crossings]]

def group_adjacent_indices(indices):
    '''
    Groups adjacent indices from the input array into separate lists.
    
    Parameters:
    indices (array-like): An array of integers representing the indices.
    
    Returns:
    list: A list of lists, where each sublist contains adjacent indices.
    
    Example:
    >>> indices = np.array([1, 2, 3, 5, 6, 8])
    >>> group_adjacent_indices(indices)
    [[1, 2, 3], [5, 6], [8]]
    '''
    # Sort the input indices to ensure they are in order
    indices = np.sort(indices)
    
    # Step 1: Find the split points based on gaps where difference > 1
    split_points = np.where(np.diff(indices) > 1)[0] + 1
    
    # Step 2: Split indices into groups based on the split points
    groups = np.split(indices, split_points)
    
    # Convert NumPy arrays in the list to Python lists
    groups = [group.tolist() for group in groups]
    
    return groups


def HeatFlowRetriveForearcMaximumCase(local_dir: str, Visit_Options):
    '''
    Retrieves the heat flow data for the forearc maximum case and saves it to a CSV file.
    
    Parameters:
    - local_dir (str): Local directory path where data is stored and outputs will be saved.
    - Visit_Options: Object containing simulation options, including all graphical times and methods to retrieve timesteps.
    
    The function:
    1. Retrieves heat flow profiles and calculates spline fits.
    2. Finds zero crossings in heat flow and identifies the forearc maximum based on the second derivative.
    3. Saves the resulting forearc maximum heat flow data to a CSV file.
    '''

    # Create directory for heat flow output if it doesn't exist
    hf_output_dir = os.path.join(local_dir, "hf_outputs")
    if not os.path.isdir(hf_output_dir):
        os.mkdir(hf_output_dir)

    # Initialize DataFrame to store results, including heat flow values and derivatives.
    df = pd.DataFrame(columns=['time', 'trench', 'longitude', 'heat flow', 'heat flux first derivative', 'heat flux second derivative'])
    
    # Loop through each graphical time in Visit_Options to retrieve heat flow profiles and process data.
    for _time in Visit_Options.all_graphical_times:
        
        # Retrieve the heat flow profile and masked longitude (Phis) for the given time.
        hfs_masked, Phis_masked = HeatFlowRetriveProfile(local_dir, _time, Visit_Options)
        phi0 = Phis_masked[0]
        phi1 = Phis_masked[-1]
        
        # Fit a smoothing spline to the heat flow data and calculate its first and second derivatives.
        hfs_spline = UnivariateSpline(Phis_masked, hfs_masked, s=0)

        # Get the corresponding timestep and vtu step for the given time.
        _time1, _, vtu_step = Visit_Options.get_timestep_by_time(_time)
        output_path = os.path.join(hf_output_dir, 'heat_flux_s%05d.csv' % vtu_step)

        # Retrieve the zero crossings and corresponding values for heat flow and its derivatives.
        zero_crossing_outputs = HeatFlowRetriveZeroCrossings(Phis_masked, hfs_spline, output_path)
        zero_crossings = zero_crossing_outputs[0]
        Phi_at_crossings = zero_crossing_outputs[1]
        hf_at_crossings = zero_crossing_outputs[2]
        dhf_dphi1_at_crossings = zero_crossing_outputs[3]
        dhf_dphi2_at_crossings = zero_crossing_outputs[4]

        # Group adjacent zero crossings indices together for processing.
        zero_crossings_groups = group_adjacent_indices(zero_crossings)
        n = 0

        # Iterate over the groups to find the forearc maximum (where second derivative < 0).
        Phi_M = np.nan; hf_M = np.nan; dhf_dphi1_M = np.nan; dhf_dphi2_M = np.nan
        for i_g in range(len(zero_crossings_groups)):
            zero_crossings_group = zero_crossings_groups[i_g]
            Phi = Phi_at_crossings[n]  # Corresponding Phi value
            dhf_dphi2 = dhf_dphi2_at_crossings[n]  # Corresponding second derivative

            print(Phi, dhf_dphi2)  # Debug print to inspect values during iteration.

            # Condition to find the forearc maximum: Phi is greater than the trench longitude and second derivative is negative.
            if Phi > (phi0 + phi1) / 2.0 and dhf_dphi2 < 0.0:
                Phi_M = Phi
                hf_M = hf_at_crossings[n]
                dhf_dphi1_M = dhf_dphi1_at_crossings[n]
                dhf_dphi2_M = dhf_dphi2_at_crossings[n]
                break
            
            n += len(zero_crossings_group)
        
        # Append the computed values for the current time step to the DataFrame.
        df = df.append({"time": _time, "trench": (phi0 + phi1) / 2.0, "longitude": Phi_M, "heat flow": hf_M,\
            "heat flux first derivative": dhf_dphi1_M, "heat flux second derivative": dhf_dphi2_M}, ignore_index=True)
    
    # Save the DataFrame to a CSV file for forearc maximum results.
    output_path = os.path.join(hf_output_dir, 'forearc_maximum.csv')
    if not os.path.isdir(os.path.basename(output_path)):
        os.mkdir(os.path.basename(output_path))
    df.to_csv(output_path, index=False)
    print(f"Output heat flux forearc maximum csv file: {output_path}")



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