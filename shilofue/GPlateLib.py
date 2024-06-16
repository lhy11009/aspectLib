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
from plate_model_manager import PlateModelManager
import pandas as pd
import gplately
import cartopy.crs as ccrs

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
    

class GPLATE_CLASS():
    '''
    A clase to work with gplate
    '''
    def __init__(self):
        '''
        initiate the class
        The **PlateModelManager** class can be used to download the plate reconstruction model files via Internet. The files will be fetched and saved into local folders. Users can use the member functions of PlateModel class to retrieve the local absolute path of the files. For example, the get_rotation_model() will return the local path to the rotation file(s).
        The **PlateReconstruction** class can be used to reconstruct tectonic plates and calculate subduction convergence stats.
        The **PlotTopologies** class can be used to plot reconstruction maps.
        Attribute:
            plate_model - a plate reconstruction model from gplate
            reconstruction_time - time of reconstruction, Ma, it must be integar
            subduction_data - reconstructed subduction zone data
            trench_ids - the valid trench ids
            age_grid_raster - reconstructed seafloor age raster
        '''
        self.anchor_plate_id = 0
        self.subduction_all_columns = ['lon', 'lat', 'conv_rate', 'conv_angle', 'trench_velocity', 'trench_velocity_angle', 'arc_length',
                                        'trench_azimuth_angle', 'subducting_pid', 'trench_pid']
        self.pm_manager = PlateModelManager()
        self.plate_model = self.pm_manager.get_model("Muller2019", data_dir="plate-model-repo")

        self.model = gplately.PlateReconstruction(self.plate_model.get_rotation_model(), 
                                             self.plate_model.get_topologies(), 
                                             self.plate_model.get_static_polygons(),
                                             anchor_plate_id=self.anchor_plate_id)
        
        self.gplot = gplately.plot.PlotTopologies(self.model, 
                                             self.plate_model.get_layer('Coastlines'), 
                                             self.plate_model.get_layer('ContinentalPolygons'), 
                                             self.plate_model.get_layer('COBs')) 

        self.reconstruction_time = 0 # initiate with reconstruction time 0

        self.subduction_data = None
        self.trench_ids = []

        self.age_grid_raster = None


    def SetReconstructionTime(self, reconstruction_time):
        '''
        set the reconstruction time
        Inputs:
            reconstruction_time - time of reconstruction
        '''
        assert(type(reconstruction_time) == int)
        self.reconstruction_time = reconstruction_time

    def Reconstruct(self):
        '''
        Reconstruct subduction zone, age raster with a reconstruction time
        # Calculate subduction convergence stats with GPlately
        # Col. 0 - longitude of sampled trench point
        # Col. 1 - latitude of sampled trench point
        # Col. 2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        # Col. 3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        # Col. 4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        # Col. 5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        # Col. 6 - length of arc segment (in degrees) that current point is on
        # Col. 7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        # Col. 8 - subducting plate ID
        # Col. 9 - trench plate ID
        '''
        # get the reconstruction of subduction zones
        self.subduction_data = self.model.tessellate_subduction_zones(self.reconstruction_time, 
                                                                # tessellation_threshold_radians=0.01, 
                                                                anchor_plate_id=self.anchor_plate_id,
                                                                ignore_warnings=True)
        # get all the trench ids
        temp = [row[9] for row in self.subduction_data]
        self.trench_ids = sorted(set(temp))

        # get the age grid raster
        self.age_grid_raster = gplately.Raster(
                                data=self.plate_model.get_raster("AgeGrids",self.reconstruction_time),
                                plate_reconstruction=self.model,
                                extent=[-180, 180, -90, 90],)
    
    def GetOneSubductionByTrenchId(self, trench_id):
        '''
        a wrapper for get_one_subduction_by_trench_id
        Inputs:
            trench_id - id of one subduction zone
        '''
        assert(trench_id in self.trench_ids)
        one_subduction_data = get_one_subduction_by_trench_id(self.subduction_data, trench_id)
        return one_subduction_data

    def PlotCoastlines(self, ax): 
        '''
        plot the coastline
        '''
        self.gplot.time = self.reconstruction_time
        self.gplot.plot_coastlines(ax, color='grey')
        

def get_one_subduction_by_trench_id(subduction_data, trench_id):
    '''
    get one subduction data from the whole dataset at one reconstruction time
    Inputs:
        subduction_data - global dataset of subduction zones at a reconstruction time
        trench_id - id of one subduction zone
    '''
    all_columns = ['lon', 'lat', 'conv_rate', 'conv_angle', 'trench_velocity', 'trench_velocity_angle', 'arc_length',
                                     'trench_azimuth_angle', 'subducting_pid', 'trench_pid']
    # ret is a list, each component is a list of 10 entries, corresponding to the 10 columns in the final panda frame 
    ret=[]
    for row in subduction_data:
        if True and (row[9]==trench_id): # only select entry with row 9
            ret.append(row)
    ret.sort(key=lambda row: row[1])
    one_subduction_data = pd.DataFrame(ret, columns=all_columns)
    
    return one_subduction_data


def ResampleSubduction(one_subduction_data, arc_length_edge, arc_length_resample_section, **kwargs):
    '''
    Resample data of a subduction to a couple of query points. 
    This is used to extract sample points of an otherwise dense subduction zone and make plots of subduction zone properties
    Inputs:
        one_subduction_data - panda data of one subduction zone
        arc_length_edge - edge of the subduction zone, filtered out
        arc_length_resample_section - resampled section of arc length
    kwargs:
        indent - indent for outputs
    '''
    # initiate
    indent = kwargs.get("indent", 0) # the default is to not indent output
    log_output_contents = ""
    data_len = len(one_subduction_data)
    
    # get the arc_length_sums
    arc_lengths = one_subduction_data['arc_length']
    arc_length_sums = np.zeros(data_len)
    arc_length_sums[0] = arc_lengths[0]
    for i in range(1, data_len):
        arc_length_sums[i] = arc_length_sums[i-1] + arc_lengths[i]

    # get the resampled arc sums
    # The resampling procedure is initiated at the center and propogated with
    # an interval of arc_length_resample_section towards both end.
    # It ends within arc_length_edge distance of the end point. 
    temp = []
    if arc_length_sums[-1] > 2 * arc_length_edge:
        temp.append(arc_length_sums[-1] / 2.0)
    i = 1
    arc_length_sum_temp = arc_length_sums[-1] / 2.0 - arc_length_resample_section
    arc_length_sum_temp1 = arc_length_sums[-1] / 2.0 + arc_length_resample_section
    while arc_length_sum_temp > arc_length_edge:
        temp.append(arc_length_sum_temp)
        temp.append(arc_length_sum_temp1)
        arc_length_sum_temp -= arc_length_resample_section
        arc_length_sum_temp1 += arc_length_resample_section
    arc_length_sums_resampled = sorted(temp)

    # resample the properties of the subduction by interpolation 
    one_subduction_data_resampled = pd.DataFrame(columns=all_columns)
    i_sbd_re = 0
    for arc_length_sum_resampled in arc_length_sums_resampled:
        for i in range(len(arc_length_sums)-1):
            if (arc_length_sums[i] <= arc_length_sum_resampled) and (arc_length_sum_resampled < arc_length_sums[i+1]):
                fraction = (arc_length_sum_resampled - arc_length_sums[i]) /  (arc_length_sums[i+1] - arc_length_sums[i])
                row_temp = fraction * one_subduction_data.iloc[i] + (1. - fraction) * one_subduction_data.iloc[i+1]
                # longitude and latitude are interpolated using a special method
                row_temp.lon, row_temp.lat = Utilities.map_mid_point(one_subduction_data.iloc[i].lon, one_subduction_data.iloc[i].lat,\
                                                                     one_subduction_data.iloc[i+1].lon, one_subduction_data.iloc[i+1].lat, fraction)
                log_output_contents += "%s%d th resampled point: (%.2f, %.2f)\n" % (" "*indent, i_sbd_re, row_temp.lon, row_temp.lat) # record the resampled point

                one_subduction_data_resampled = pd.concat([one_subduction_data_resampled,  pd.DataFrame([row_temp])], ignore_index=True)
        i_sbd_re += 1

    return one_subduction_data_resampled, log_output_contents


def plot_one_subduction_data(ax, one_subduction_data, **kwargs):
    '''
    plot one subduction zone in a global projection
    Inputs:
        ax - plot axis
        one_subduction_data - dataset of one subduction zone
    '''
    attribute = kwargs.get("attribute", "conv_rate") # the attribute to plot
    # Latitudes and longitudes of points along trench segments
    subduction_lon = one_subduction_data.lon
    subduction_lat = one_subduction_data.lat

    # get the data to plot with
    data_attr = getattr(one_subduction_data, attribute) 

    # figure configurations
    vmin = None; vmax = None
    if attribute == "conv_rate":
        vmin = 0.0
        vmax = 20.0
        cmap="inferno_r"
    elif attribute == 'trench_velocity':
        vmin = -5.0
        vmax = 5.0
        cmap="viridis_r"
    else:
        raise NotImplementedError()
    
    # plot
    cb=ax.scatter(subduction_lon,subduction_lat, marker=".", s=5, c=data_attr, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)

    return cb


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