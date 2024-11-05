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
from multiprocessing.sharedctypes import Value
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
import re
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
    
    def GetSubductionData(self):
        '''
        Return subduction data as pandas frame
        '''
        return pd.DataFrame(self.subduction_data, columns=self.subduction_all_columns)
    
    def GetOneSubductionByTrenchId(self, trench_id):
        '''
        a wrapper for get_one_subduction_by_trench_id
        Inputs:
            trench_id - id of one subduction zone
        '''
        assert(trench_id in self.trench_ids)
        one_subduction_data = get_one_subduction_by_trench_id(self.subduction_data, trench_id)
        return one_subduction_data
    
    def ResampleSubductionById(self, trench_id, arc_length_edge, arc_length_resample_section):
        '''
        Resample an existing subduction zone by id
        Inputs:
            trench_id - id of one subduction zone
        '''
        assert(trench_id in self.trench_ids)
        one_subduction_data = get_one_subduction_by_trench_id(self.subduction_data, trench_id)
        one_subduction_data_resampled, log_output_contents = resample_subduction(one_subduction_data, arc_length_edge, arc_length_resample_section, self.subduction_all_columns, indent=4)
        return one_subduction_data_resampled, log_output_contents
    
    def ResampleAllSubduction(self, arc_length_edge, arc_length_resample_section):
        '''
        Resample an existing subduction zone by id
        '''
        # subduction_data_resampled = pd.DataFrame(columns=self.subduction_all_columns)
        subduction_data_resampled = None
        log_output_contents = ""
        for i in range(len(self.trench_ids)):
            trench_id = self.trench_ids[i]
            one_subduction_data = get_one_subduction_by_trench_id(self.subduction_data, trench_id)
            one_subduction_data_resampled, log_output_contents = resample_subduction(one_subduction_data, arc_length_edge, arc_length_resample_section, self.subduction_all_columns, indent=4)
            log_output_contents += "%d th arc\n" % i
            log_output_contents += "start (%.2f, %.2f)\n" % (one_subduction_data.iloc[0].lon, one_subduction_data.iloc[0].lat)
            log_output_contents += "end (%.2f, %.2f)\n" % (one_subduction_data.iloc[-1].lon, one_subduction_data.iloc[-1].lat)
            # log_output_contents += "length in degree: %.2f\n" % (one_subduction_data["arc_length"][len(one_subduction_data)-1])
            if i == 0:
                # subduction_data_resampled = pd.DataFrame([one_subduction_data_resampled])
                subduction_data_resampled = pd.DataFrame(one_subduction_data_resampled)
            else:
                # subduction_data_resampled = pd.concat([subduction_data_resampled,  pd.DataFrame([one_subduction_data_resampled])], ignore_index=True)
                subduction_data_resampled = pd.concat([subduction_data_resampled,  pd.DataFrame(one_subduction_data_resampled)], ignore_index=True)
        return subduction_data_resampled

    def InterpolateAgeGrid(self, samples):
        '''
        InterpolateAgeGrid at sample points
        ''' 
        ages = self.age_grid_raster.interpolate(samples.lon, samples.lat, method="nearest")
        return ages

    def PlotCoastlines(self, ax): 
        '''
        plot the coastline
        '''
        self.gplot.time = self.reconstruction_time
        self.gplot.plot_coastlines(ax, color='grey')

    def PlotSeaFloorAges(self, ax):
        '''
        plot the seafloor ages
        '''
        im_age = self.gplot.plot_grid(ax, self.age_grid_raster.data, cmap='YlGnBu', vmin=0, vmax=200, alpha=0.8)
        return im_age

    def FixTrenchAge(self, subduction_data, **kwargs):
        '''
        Fix the trench ages in subduction_data
        Inputs:
            subduction_data: pandas object, subduction dataset
        '''
        # automatically fix the invalid ages 
        for i in range(len(subduction_data)):
            fix_age_polarity = subduction_data.fix_age_polarity[i]
            if not np.isnan(fix_age_polarity):
                # fix with existing polarity
                # 0 and 1: on different side of the trench
                # 2: manually assign values of longitude and latitude
                if (fix_age_polarity == 0): 
                    new_age = self.FixTrenchAgeLocal(subduction_data, i, subduction_data.trench_azimuth_angle[i] + 180.0)
                elif (fix_age_polarity == 1): 
                    new_age = self.FixTrenchAgeLocal(subduction_data, i, subduction_data.trench_azimuth_angle[i])
                elif (fix_age_polarity == 2):
                    subduction_data_local0 = pd.DataFrame([subduction_data.iloc[i]])
                    subduction_data_local0.lon, subduction_data_local0.lat = subduction_data.iloc[i].lon_fix, subduction_data.iloc[i].lat_fix
                    new_age = self.InterpolateAgeGrid(subduction_data_local0)
                    subduction_data['age'][i] = new_age
                    pass
                else:
                    raise NotImplementedError
            else:
                # figure out a possible polarity
                new_age = self.FixTrenchAgeLocal(subduction_data, i, subduction_data.trench_azimuth_angle[i] + 180.0)
                if np.isnan(new_age):
                    # next, try the other direction
                    new_age = self.FixTrenchAgeLocal(subduction_data, i, subduction_data.trench_azimuth_angle[i])
                    if not np.isnan(new_age):
                        subduction_data["fix_age_polarity"][i] = 1
                else:
                    subduction_data["fix_age_polarity"][i] = 0

    
    def FixTrenchAgeLocal(self, subduction_data, i, theta):
        '''
        fix the invalid age in a subduction data object with the age grid in the class
        Inputs:
            subduction_data - subduction data
            i - ith index to fix
            theta - direction to look for new point
        '''
        ds = [12.5e3, 25e3, 50e3, 75e3, 100e3, 150e3, 200e3, 300e3, 400e3]
        new_age = np.nan
        for j in range(len(ds)-1):
            # first get two points
            subduction_data_local0 = pd.DataFrame([subduction_data.iloc[i]])
            subduction_data_local1 = pd.DataFrame([subduction_data.iloc[i]])
            subduction_data_local0.lon, subduction_data_local0.lat = \
                Utilities.map_point_by_distance(subduction_data.iloc[i].lon, subduction_data.iloc[i].lat, theta, ds[j])
            subduction_data_local1.lon, subduction_data_local1.lat = \
                Utilities.map_point_by_distance(subduction_data.iloc[i].lon, subduction_data.iloc[i].lat, theta, ds[j+1])
            # then interpolate ages at these two points
            new_age0 = self.InterpolateAgeGrid(subduction_data_local0)
            new_age1 = self.InterpolateAgeGrid(subduction_data_local1)
            # if both points are valid, intepolate with these two points
            if (not np.isnan(new_age0)) and (not np.isnan(new_age1)):
                # new_age = (new_age0 * (x - ds[j+1]) + new_age1 *(x - ds[j])) / (ds[j+1] - ds[j]) and x = 0.0
                new_age = (new_age0 * ds[j+1] - new_age1 * ds[j]) / (ds[j+1] - ds[j])
                # print("i = %d, new_age = %.2f" % (i, new_age)) # debug
                subduction_data.age.iloc[i] = new_age
                subduction_data.lon_fix.iloc[i] = subduction_data_local1.lon # note this records the further point in this case
                subduction_data.lat_fix.iloc[i] = subduction_data_local0.lat
                break
            else:
                subduction_data.age.iloc[i] = np.nan
        return new_age


class PARSERECONSTRUCTION():
    def __init__(self):
        '''
        initiation
        '''
        self.trench_data = []  # record the trench coordinates
        self.trench_names = []
        self.trench_pids = []
        self.trench_begin_times = []
        self.trench_end_times = []
        pass

    def ReadFile(self, infile):
        '''
        read a file of reconstruction data
        Inputs:
            infile - input file path
        '''
        infile = os.path.join(ASPECT_LAB_DIR, "dtemp", "gplate_export_test0", "Muller_etal_2019_PlateBoundaries_no_topologies", "reconstructed_0.00Ma.xy")
        assert(os.path.isfile(infile))
        # Read a file of plate reconstruction.
        # We'll look for the key word "SubductionZone" in the file
        # and save there respectvie coordinates into a big array
        i=0
        temp_l = []
        temp_d = []
        n_trench = 0  # record the number of trenches
        sbd_begin = False  # tags for in a section of trench data
        sbd_end = False
        read=True # a tag for continuing reading the file
        with open(infile, 'r') as fin:
            line = fin.readline()
            i += 1
            while line:
                read = True # by default, we will continue to read the file in each loop
                # First, in case we are reading a data section of trench locations,
                # check if we already reached the end of a section
                if sbd_begin and re.match('^>', line):
                    sbd_end = True
                # Then, determine what to do.
                # 1. we are reading a data section and want to continue
                # 2. we are reading a data section and want to end
                # 3. we are not reading a data section and is haunting for the next
                # section in the file.
                if sbd_begin and (not sbd_end):
                    # begin read in subduction zone data
                    # To convert string inputs to data:
                    #   first, split the string by any number of whitespace characters
                    #   then, convert the data to floats and all the additional whitespace
                    #   shoule be removed
                    temp_data = line.split()
                    temp_data = [float(x) for x in temp_data]
                    temp_d.append(temp_data)
                    pass
                elif sbd_begin and sbd_end:
                    # reach the end of a section
                    self.trench_data.append(temp_d)
                    sbd_begin = False  # reset the two tags
                    sbd_end = False
                    read = False
                    pass
                elif re.match('^>SubductionZone', line):
                    # these inputs at the begining of a line
                    # indicate the start of a subduction zone
                    # data section.
                    temp_l.append(i)
                    sbd_begin = True # set the tag for the begining of the dataset
                    temp_d = []
                    # after finding the start of the section,
                    # continue to read through the headers
                    # of this section
                    while(line and re.match('^>', line)):
                        line=fin.readline()
                        i += 1
                        if re.match('^> name', line):
                            # look for the name of the trench and append it to the list
                            self.trench_names.append(Utilities.remove_substrings(line, ["> name ", '\n']))
                        elif re.match('> reconstructionPlateId', line):
                            # look for the plate id of the trench and append it to the list
                            self.trench_pids.append(int(Utilities.remove_substrings(line, ["> reconstructionPlateId ", '\n'])))
                        elif re.match('> validTime TimePeriod <begin> TimeInstant <timePosition>', line):
                            # look for the begin and end time of the trench and append it to the list
                            temp0 = Utilities.remove_substrings(line, ["> validTime TimePeriod <begin> TimeInstant <timePosition>", '</timePosition>.*\n'])
                            self.trench_begin_times.append(float(temp0))
                            temp1 = Utilities.remove_substrings(line, ['^.*<end> TimeInstant <timePosition>', '</timePosition>.*\n'])
                            if type(temp1) == float:
                                self.trench_end_times.append(temp1)
                            else:
                                self.trench_end_times.append(0.0)
                    read = False
                if read:
                    # if the flag marks continue reading the file, proceed
                    line = fin.readline()
                    i += 1
        i -= 1 # the last time reading is unsuccessful
        n_trench = len(self.trench_data)

    def LookupNameByPid(self, pid):
        '''
        Inputs:
            pid - index of the plate
        '''
        _name = ""
        assert(type(pid) == int)
        try:
            _index = self.trench_pids.index(pid)
        except ValueError:
            _name = ""
        else:
            _name = self.trench_names[_index]
        return _name
    
        
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


def resample_subduction(one_subduction_data, arc_length_edge, arc_length_resample_section, all_columns, **kwargs):
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
    ssize = kwargs.get('s', 5)
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
    im = ax.scatter(subduction_lon,subduction_lat, marker=".", s=ssize, c=data_attr, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)

    return im


def MaskBySubductionTrenchIds(subduction_data, subducting_pid, trench_pid, i_p):
    """
    Generates a combined mask for subduction data based on user selection or specific 
    subducting and trench IDs.
    
    Parameters:
        subduction_data (pd.DataFrame): The DataFrame containing subduction data to be filtered.
        subducting_pid (int or None): The subducting plate ID to match. If None, all IDs are included.
        trench_pid (int or None): The trench plate ID to match. If None, all IDs are included.
        i_p (list or None): List of indices selected by the user. If not None, these indices are used.
    
    Returns:
        np.ndarray: A boolean mask combining the specified conditions for filtering the data.
    
    Implementation:
        - If `i_p` is provided, create `mask1` to select only those indices.
        - If `subducting_pid` is provided, create `mask1` to select rows matching the `subducting_pid`.
        - If neither is provided, `mask1` includes all rows.
        - If `trench_pid` is provided, create `mask2` to select rows matching the `trench_pid`.
        - If `trench_pid` is not provided, `mask2` includes all rows.
        - The final mask is the logical AND of `mask1` and `mask2`.
    """
    if i_p is not None:
        mask1 = np.zeros(len(subduction_data), dtype=bool)
        mask1[i_p] = 1
    elif subducting_pid is not None:
        # Generate mask1 based on the provided subducting plate ID
        mask1 = subduction_data.subducting_pid == subducting_pid
    else:
        mask1 = np.ones(len(subduction_data), dtype=bool)

    if trench_pid is not None:
        # Generate mask2 based on the provided trench plate ID
        mask2 = subduction_data.trench_pid == trench_pid
    else:
        mask2 = np.ones(len(subduction_data), dtype=bool)

    return (mask1 & mask2)


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