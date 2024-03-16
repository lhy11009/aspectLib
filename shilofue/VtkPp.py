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
# from matplotlib import cm
from matplotlib import pyplot as plt
import vtk
import scipy.interpolate
from vtk.util.numpy_support import vtk_to_numpy
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities
import time

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


class VERTICAL_PROFILE():
    '''
    A vertical profile
    Attributes:
        rs - vertical coordinates
        n - points in the profile
        field_funcs (dict) - functions in the profile
        uniform - is the grid uniform
    '''
    def __init__(self, rs, field_names, field_datas, uniform=False):
        '''
        Initiation
        '''
        self.rs = rs
        self.n = rs.size
        self.uniform = uniform
        self.field_funcs = {}
        for i in range(len(field_names)):
            field_name = field_names[i]
            field_data = field_datas[i]
            field_func = scipy.interpolate.interp1d(rs, field_data, assume_sorted=True, fill_value="extrapolate")
            self.field_funcs[field_name] = field_func
        pass

    def GetFunction(self, field_name):
        '''
        field_name(str): name of the field
        return
            function of the matched field
        '''
        return self.field_funcs[field_name]


class VTKP():
    '''
    Class for vtk post-process utilities
    Attributes:
        reader (vtk reader object): reader for vtk file format
        i_poly_data (vtkPolydata): input data
    '''
    def __init__(self, **kwargs):
        '''
        Initiation
        Inputs:
            kwargs:
                dim - dimension
        '''
        self.i_poly_data = vtk.vtkPolyData()
        self.include_cell_center = False
        self.c_poly_data = vtk.vtkPolyData()
        self.cell_sizes = None
        self.dim = kwargs.get('dim', 2)
        self.grav_data = None  # a 2 column array to save the gravity data (depth in meter and grav_acc)
        self.geometry = kwargs.get('geometry', 'chunk')
        self.Ro = kwargs.get('Ro', 6371e3)
        self.Xmax = kwargs.get('Xmax', 61.0 * np.pi / 180.0)
        self.time = kwargs.get('time', None)
        ha_file = kwargs.get("ha_file", None)
        if ha_file is None:
            self.Tref_func = None
        else:
            assert(os.path.isfile(ha_file))
            DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
            DepthAverage.Import(ha_file)
            Utilities.my_assert(self.time != None, ValueError, "\"time\" is a requried input if \"ha_file\" is presented")
            self.Tref_func = DepthAverage.GetInterpolateFunc(self.time, "temperature")
            self.density_ref_func = DepthAverage.GetInterpolateFunc(self.time, "adiabatic_density")
        pass

    def ReadFile(self, filein):
        '''
        Read file
        Inputs:
            filein (str): path to a input file
        '''
        assert(os.path.isfile(filein))
        file_extension = filein.split('.')[-1]
        if file_extension == 'pvtu':
            self.reader = vtk.vtkXMLPUnstructuredGridReader()
        else:
            raise TypeError("%s: Wrong type of file" % Utilities.func_name())
        self.reader.SetFileName(filein)
        self.reader.Update()
        print("ReadFile: %s" % filein)
    
    def ConstructPolyData(self, field_names, **kwargs):
        '''
        construct poly data
        Inputs:
            field_names (list of field names)
            kwargs:
                include_cell_center - include_cell_center in the poly_data
        '''
        start = time.time()
        include_cell_center = kwargs.get('include_cell_center', False)
        fix_cell_value = kwargs.get('fix_cell_value', True)  # fix the value of cell centers
        quiet = kwargs.get('quiet', False)
        construct_Tdiff = kwargs.get("construct_Tdiff", False)
        assert(type(field_names) == list and len(field_names) > 0)
        noP = 0
        noC = 0
        grid = self.reader.GetOutput()
        data_set = self.reader.GetOutputAsDataSet()
        points = grid.GetPoints()
        cells = grid.GetCells()
        point_data = data_set.GetPointData()
        self.i_poly_data.SetPoints(points)
        self.i_poly_data.SetPolys(cells)
        # import polydata
        is_first = True
        for field_name in field_names:
            if is_first:
                self.i_poly_data.GetPointData().SetScalars(point_data.GetArray(field_name))  # put T into cell data
                is_first = False
            else:
                self.i_poly_data.GetPointData().AddArray(point_data.GetArray(field_name))  # put T into cell data
            noP = self.i_poly_data.GetNumberOfPoints()
        if construct_Tdiff:
            assert(self.Tref_func != None)
            self.ConstructTemperatureDifference()
        time_import = time.time()
        if include_cell_center:
            centers = vtk.vtkCellCenters()  # test the utilities of centers
            centers.SetInputData(self.i_poly_data)
            centers.Update()
            probeFilter = vtk.vtkProbeFilter()
            # probeFilter = vtk.vtkCompositeDataProbeFilter() # debug
            probeFilter.SetSourceData(self.i_poly_data)  # use the polydata
            probeFilter.SetInputData(centers.GetOutput()) # Interpolate 'Source' at these points
            probeFilter.Update()
            self.c_poly_data = probeFilter.GetOutput()  # poly data at the center of the point
            self.include_cell_center = True
            noC = self.c_poly_data.GetNumberOfPoints() 
            # self.c_poly_data.GetPointData().GetArray('T') = numpy_to_vtk(T_field)
            # cell sizes
            cell_size_filter = vtk.vtkCellSizeFilter()
            cell_size_filter.SetInputDataObject(grid)
            cell_size_filter.Update()
            cz_cell_data = cell_size_filter.GetOutputDataObject(0).GetCellData()
            if self.dim == 2:
                self.cell_sizes = vtk_to_numpy(cz_cell_data.GetArray('Area'))
            else:
                raise ValueError("Not implemented")
            # fix values of fields. This is needed because the interpolation is not correct where 
            # the mesh refines or coarsens. 
            # I am usign the 'T' field as an indicator:
            # every cell center with T = 0.0 is to be fixed (assuming a realistic T > 273.0)
            # the strategy is to take a nearby cell center and check its value.
            # Note: this will look into the adjacent cells until it finds one with a sufficently 
            # approximate location and a non-zero value of T.
            tolerance = 1.0
            T_field = vtk_to_numpy(self.c_poly_data.GetPointData().GetArray('T'))
            fields = []
            for field_name in field_names:
                fields.append(vtk_to_numpy(self.c_poly_data.GetPointData().GetArray(field_name)))
            # density_field =  vtk_to_numpy(self.c_poly_data_raw.GetPointData().GetArray('density'))
            if fix_cell_value:
                for i in range(noC):
                    if T_field[i] - 0.0 < tolerance:
                        xs = self.c_poly_data.GetPoint(i)
                        found = False
                        i1 = 0
                        j = 1
                        dist_max = 3*(self.cell_sizes[i]**0.5)  # compare to the cell size
                        while True:   # find a adjacent point
                            if i+j >= noC and i-j < 0:
                                break # end is reached
                            if i+j < noC:
                                xs1 = self.c_poly_data.GetPoint(i+j)
                                dist = ((xs1[0] - xs[0])**2.0 + (xs1[1] - xs[1])**2.0)**0.5
                                if T_field[i+j] - 0.0 > tolerance and dist < dist_max:
                                    i1 = i + j
                                    found = True
                                    break
                            if i-j >= 0:
                                xs1 = self.c_poly_data.GetPoint(i-j)
                                dist = ((xs1[0] - xs[0])**2.0 + (xs1[1] - xs[1])**2.0)**0.5
                                if i-j >= 0 and T_field[i-j] - 0.0 > tolerance and dist < dist_max:
                                    i1 = i - j
                                    found = True
                                    break
                            j += 1
                        if not found:
                            raise ValueError("A cell center (%.4e, %.4e, %.4e) is not in mesh, and the adjacent cells are not found" % (xs[0], xs[1], xs[2]))
                        for n in range(len(fields)):
                            fields[n][i] = fields[n][i1]
        time_center = time.time()
        # send out message
        message = "ConstructPolyData: %d * (%d + %d) entries in the polydata imported and %d * (%d + %d) points in the data at cell center. \
import data takes %f, interpolate cell center data takes %f"\
        % (noP, self.dim, len(field_names), noC, self.dim, len(field_names), time_import - start, time_center - time_import)
        if not quiet:
            print(message)


    def ConstructTemperatureDifference(self):
        '''
        Construct a dT field of temperature differences
        '''
        T_field = vtk_to_numpy(self.i_poly_data.GetPointData().GetArray("T"))
        Tdiffs = vtk.vtkFloatArray()
        Tdiffs.SetName("dT")
        for i in range(self.i_poly_data.GetNumberOfPoints()):
            xs = self.i_poly_data.GetPoint(i)
            x =  xs[0]
            y =  xs[1]
            r = get_r(x, y, self.geometry)
            Tref = self.Tref_func(self.Ro - r)
            T = T_field[i]
            Tdiffs.InsertNextValue(T - Tref)
        self.i_poly_data.GetPointData().AddArray(Tdiffs)

    def VerticalProfile2D(self, x0_range, x1, n, **kwargs):
        '''
        Get vertical profile by looking at one vertical line
        Inputs:
            x0_range: range of the first coordinate (y or r)
            x1: second coordinate (x or theta)
            n (int): number of point in the profile
            kwargs (dict):
                geometry: spherical or cartesian
                fix_point_value: fix invalid value indicated by zero temperature
                                these points exist as a result of interpolation
                                in adaptive meshes.
        '''
        geometry = kwargs.get('geometry', 'chunk')
        fix_point_value = kwargs.get('fix_point_value', False)
        assert(len(x0_range) == 2)
        points = np.zeros((n, 2))
        x0s = np.linspace(x0_range[0], x0_range[1], n)
        x1s = np.ones(n) * x1
        if geometry == 'chunk':
            xs = x0s * np.cos(x1s)
            ys = x0s * np.sin(x1s)
        elif geometry == 'box':
            ys = x0s
            xs = x1s
        else:
            raise ValueError('Wrong option for geometry')
        points[:, 0] = xs
        points[:, 1] = ys
        v_poly_data = InterpolateGrid(self.i_poly_data, points, quiet=True)
        point_data = v_poly_data.GetPointData()
        fields = []
        field_names = []
        fix_value_mask = None
        for i in range(point_data.GetNumberOfArrays()):
            field_name = point_data.GetArrayName(i)
            field = vtk_to_numpy(point_data.GetArray(field_name))
            field_names.append(field_name)
            fields.append(field)
            if fix_point_value and field_name == "T":
                # only take the valid point value (T = 0.0)
                fix_value_mask = (field > 1e-6)
        fields1 = []
        if fix_value_mask is not None:
            x0s = x0s[fix_value_mask]
            for field in fields:
                fields1.append(field[fix_value_mask])
        else:
            fields1 = fields
        v_profile = VERTICAL_PROFILE(x0s, field_names, fields1, uniform=True)
        return v_profile
    
    def StaticPressure(self, x0_range, x1, n, **kwargs):
        '''
        Compute the static pressure at a point
        Inputs:
            x0_range: range of the first coordinate (y or r)
            x1: second coordinate (x or theta)
            n (int): number of point in the profile
            kwargs (dict):
                geometry: spherical or cartesian
        '''
        geometry = kwargs.get('geometry', 'chunk') # geometry
        use_gravity_profile = False
        if self.grav_data is not None:
            use_gravity_profile = True
        else:
            constant_grav_acc = kwargs['grav_acc']
        points = np.zeros((n, 2))
        assert(x0_range[0] < x0_range[1])
        x0s = np.linspace(x0_range[0], x0_range[1], n)
        interval = (x0_range[1] - x0_range[0]) / (n-1)
        x1s = np.ones(n) * x1
        if geometry == 'chunk':
            xs = x0s * np.cos(x1s)
            ys = x0s * np.sin(x1s)
        elif geometry == 'box':
            ys = x0s
            xs = x1s
        else:
            raise ValueError('Wrong option for geometry')
        points[:, 0] = xs
        points[:, 1] = ys
        v_poly_data = InterpolateGrid(self.i_poly_data, points, quiet=True)
        point_data = v_poly_data.GetPointData()
        density_field = vtk_to_numpy(point_data.GetArray('density'))
        static_pressure = 0.0  # integrate the static pressure
        for i in range(1, n):
            depth = self.Ro - (x0s[i] + x0s[i-1])/2.0
            if use_gravity_profile:
                grav_acc = self.GetGravityAcc(depth)
            else:
                grav_acc = constant_grav_acc # geometry
            static_pressure += (density_field[i-1] + density_field[i]) / 2.0 * grav_acc * interval
        return static_pressure
    
    def ImportGravityData(self, filein):
        '''
        Import gravity data, file shoud contain depth and 
        gravity accerleration.
        Inputs:
            filein (str): path to a input file
        '''
        assert(os.path.isfile(filein))
        self.grav_data = np.loadtxt(filein, skiprows=8)
    
    def GetGravityAcc(self, depths):
        '''
        Get gravity from a profile
        Inputs:
            depths - depth of point, float or a vector
        '''
        assert(self.grav_data.shape[1] == 2)
        grav_accs = np.interp(depths, self.grav_data[:, 0], self.grav_data[:, 1])
        return grav_accs


# todo_export
def ExportPointGridFromPolyData(i_poly_data, ids, output_xy=False):
    '''
    export point grid from a given vtk poly data by the indexed
    '''
    assert(ids is not [])
    vtk_points = vtk.vtkPoints()
    for id in ids:
        xs = i_poly_data.GetPoint(id)
        vtk_points.InsertNextPoint(xs[0], xs[1], xs[2])
    point_grid = vtk.vtkUnstructuredGrid()
    point_grid.SetPoints(vtk_points)
    if output_xy:
        coords = vtk_to_numpy(point_grid.GetPoints().GetData())
        return coords
    else:
        return point_grid
    pass


def ExportContour(poly_data, field_name, contour_value, **kwargs):
    '''
    Export contour of a field with a value
    Inputs:
        field_name (str): name of the field
        contour_value (float): contour value
        kwargs:
            fileout (str): path to output
    '''
    print("Filter contour")
    fileout = kwargs.get('fileout', None)
    contour_filter = vtk.vtkContourFilter()
    # prepare poly data for contour
    c_poly_data = vtk.vtkPolyData()
    c_vtk_point_data = poly_data.GetPointData()  # vtkPointData
    c_poly_data.SetPoints(poly_data.GetPoints())  # import points and polys
    c_poly_data.SetPolys(poly_data.GetPolys())
    vtk_data_array = c_vtk_point_data.GetArray(field_name)
    assert(vtk_data_array != None)
    c_poly_data.GetPointData().SetScalars(vtk_data_array)
    # draw contour 
    contour_filter.SetInputData(c_poly_data)
    contour_filter.Update()
    contour_filter.GenerateValues(1, contour_value, contour_value)  # Extract just one contour
    contour_filter.Update()
    # write output if a path is provided
    if fileout != None:
        file_extension = fileout.split('.')[-1]
        if file_extension == 'txt':
            writer = vtk.vtkSimplePointsWriter()
        else:
            raise TypeError("%s: Wrong type of file" % Utilities.func_name())
        writer.SetInputData(contour_filter.GetOutput())
        writer.SetFileName(fileout)
        writer.Update()
        writer.Write()
        print("%s, Generate output (contour) : %s" % (Utilities.func_name(), fileout))
    return contour_filter.GetOutput()


def ExportPolyData(poly_data, fileout, **kwargs):
    '''
    Export poly data to vtp file
    '''
    indent = kwargs.get('indent', 0)
    # output
    file_extension = fileout.split('.')[-1]
    if file_extension == "vtp":
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileout)
        writer.SetInputData(poly_data)
        # writer.SetFileTypeToBinary()  # try this later to see if this works
        writer.Update()
        writer.Write()
    elif file_extension == 'txt':
        raise TypeError("%s: option for txt file is not implemented" % Utilities.func_name())
        # ExportPolyDataAscii(poly_data, fileout)
    else:
        raise TypeError("%s: Wrong type of file" % Utilities.func_name())
    print(' '*indent + "%s: Write file %s" % (Utilities.func_name(), fileout))


def ExportPolyDataAscii(poly_data, field_names, file_out):
    '''
    Export Poly data to a ascii file
    '''
    print("%s: operating" % Utilities.func_name())
    header = "# 1: x (m)\n# 2: y (m)"
    i = 3
    point_data_export = []  # point data
    for field_name in field_names:
        header += "\n# %d" % i + ": " + field_name
        vtk_data_array = poly_data.GetPointData().GetArray(field_name)
        Utilities.my_assert(vtk_data_array != None, KeyError,\
            "Failed to get field %s from vtk point data" % field_name)
        point_data_export.append(vtk_to_numpy(vtk_data_array))
        i += 1
    header += '\n'
    # output data
    output= ""
    for i in range(poly_data.GetNumberOfPoints()):
        if i > 0:
            output += "\n"  # append new line for new point
        xs = poly_data.GetPoint(i)
        output += "%.4e" % xs[0] + "\t" + "%.4e" % xs[1]
        j = 0
        for field_name in field_names:
            val = point_data_export[j][i]; 
            output += "\t" + "%.4e" % val
            j += 1
    # write file
    with open(file_out, 'w') as fout:
        fout.write(header)
        fout.write(output)
    print("\tWrite ascii data file: %s" % file_out)


def ExportPolyDataFromRaw(Xs, Ys, Zs, Fs, fileout, **kwargs):
    '''
    Export poly data from raw data
    '''
    i_points = vtk.vtkPoints()
    field_name = kwargs.get("field_name", "foo")
    assert(Xs.size == Ys.size)
    if Zs != None:
        assert(Xs.size == Zs.size)
    for i in range(Xs.size):
        x = Xs[i]
        y = Ys[i]
        if Zs is not None:
            z = Zs[i]
        else:
            z = 0.0
        i_points.InsertNextPoint(x, y, z)
    i_poly_data = vtk.vtkPolyData()  # initiate poly daa
    i_poly_data.SetPoints(i_points) # insert points
    if Fs != None:
        # insert field data
        assert(Xs.size == Fs.size)
        fvalues.SetName(field_name)
        i_poly_data.GetPointData().SetScalars(numpy_to_vtk(Fs, deep=1))
    ExportPolyData(i_poly_data, fileout, **kwargs)


def InterpolateGrid(poly_data, points, **kwargs):
    '''
    Inputs:
        poly_data (vtkPolyData): input data set
        points (np array): grid to interpolate to
    Return:
        poly data on the new grid
    Output:
    '''
    quiet = kwargs.get('quiet', False)
    fileout = kwargs.get('fileout', None)
    assert(points.ndim == 2)  # points is a 2d array
    if not quiet:
        print("%s: Perform interpolation onto the new grid" % Utilities.func_name())
    grid_points = vtk.vtkPoints()
    for i in range(points.shape[0]):
        x = points[i][0]
        y = points[i][1]
        grid_points.InsertNextPoint(x, y, 0)  # this is always x, y, z
    grid_data = vtk.vtkPolyData()
    grid_data.SetPoints(grid_points)
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(poly_data)  # use the polydata
    probeFilter.SetInputData(grid_data) # Interpolate 'Source' at these points
    probeFilter.Update()
    # export to file output if a valid file is provided
    if fileout != None:
        print("\t%s: export data" % Utilities.func_name())
        ExportPolyData(probeFilter.GetOutput(), fileout, indent=4)
    return probeFilter.GetOutput()


def NpIntToIdList(numpy_int_array):
    '''
    Convert numpy in array to vtk id list
    Inputs:
        numpy_int_array - 1d numpy int array
    '''
    id_list = vtk.vtkIdList()
    if type(numpy_int_array) == list:
        for i in range(len(numpy_int_array)):
            id_list.InsertNextId(numpy_int_array[i])
    else:
        assert(numpy_int_array.ndim == 1)  # 1d array
        for i in range(numpy_int_array.size):
            id_list.InsertNextId(numpy_int_array[i])
    return id_list


def OperateDataArrays(poly_data, names, operations):
    '''
    perform operation to vtk_arrays
    Inputs:
        poly_data - a poly data to work with
        names - names of fields to operate on
        operations - 0 (add) or 1 (minus)
    '''
    assert(len(operations) == len(names) - 1)
    is_first = True
    i = 0
    for _name in names:
        if is_first:
            o_array = np.copy(vtk_to_numpy(poly_data.GetArray(_name)))
            is_first = False
        else:
            if operations[i] == 0:
                o_array += vtk_to_numpy(poly_data.GetArray(_name))
            elif operations[i] == 1:
                o_array -= vtk_to_numpy(poly_data.GetArray(_name))
            else:
                raise ValueError('operation for %d is not implemented.' % operations[i])
            i += 1
    return o_array


def get_r(x, y, geometry): 
    '''
    Get r (the first coordinate)
    Inputs:
        x - x coordinate
        y - y coordinate
        geometry - 'chunk' or 'box'
    '''
    if geometry == 'chunk':
        r = (x*x + y*y)**0.5
    elif geometry == 'box':
        r = y
    else:
        raise ValueError("not implemented")
    return r


def get_r3(x, y, z, geometry):
    '''
    Get r (the first coordinate)
    Inputs:
        x - x coordinate
        y - y coordinate
        geometry - 'chunk' or 'box'
    '''
    if geometry == 'chunk':
        r = (x*x + y*y + z*z)**0.5
    elif geometry == 'box':
        r = z
    else:
        raise ValueError("not implemented")
    return r


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