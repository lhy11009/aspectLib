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
from vtk.util.numpy_support import vtk_to_numpy

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

class VTKP():
    '''
    Class for vtk post-process utilities
    Attributes:
        reader (vtk reader object): reader for vtk file format
        i_poly_data (vtkPolydata): input data
    '''

    def __init__(self):
        '''
        Initiation
        '''
        reader = vtk.vtkXMLPUnstructuredGridReader()
        self.i_poly_data = vtk.vtkPolyData()
        pass

    def ReadFile(self, filein):
        '''
        Read file
        Inputs:
            filein (str): path to a input file
        '''
        filein = "./solution/solution-00002.pvtu"
        assert(os.path.isfile(filein))
        self.reader.SetFileName(filein)
        self.reader.Update()
    
    def ConstructPolyData(self, field_names):
        '''
        construct poly data
        Inputs:
            field_names (list of field names)
        '''
        assert(type(field_names) == list and len(field_names) > 0)
        grid = self.reader.GetOutput()
        data_set = self.eader.GetOutputAsDataSet()
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
    c_poly_data.SetPoints(poly_data.GetPoints())
    c_poly_data.GetPointData().SetScalars(c_vtk_point_data.GetArray(field_name))
    # draw contour 
    contour_filter.SetInputData(c_poly_data)
    contour_filter.Update()
    contour_filter.GenerateValues(1, contour_value, contour_value)  # Extract just one contour
    contour_filter.Update()
    # write output if a path is provided
    if (os.access(fileout, os.R_OK)):
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


def ExportPolyData(poly_data, fileout):
    '''
    Export poly data to vtp file
    '''
    assert(os.access(fileout, os.R_OK))
    # output
    file_extension = fileout.split('.')[-1]
    if file_extension == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    else:
        raise TypeError("%s: Wrong type of file" % Utilities.func_name())
    writer.SetFileName(fileout)
    writer.SetInputData(poly_data)
    # writer.SetFileTypeToBinary()  # try this later to see if this works
    writer.Update()
    writer.Write()
    print("%s: Write file %s" % (Utilities.func_name(), fileout))


def ExportPolyDataAscii(poly_data, fileout, field_names):
    '''
    Export Poly data to a ascii file
    '''
    print("%s: operating" % Utilities.func_name())
    header = "# 1: x (m)\n# 2: y (m)"
    i = 3
    for field_name in field_names:
        header += "\n# %d" % i + ": " + field_name
        i += 1
    # output data
    output= ""
    for i in range(poly_data.GetNumberOfPoints()):
        if i > 0:
            output += "\n"  # append new line for new point
        xs = poly_data.GetPoint(i)
        output += xs[0] + "\t" + xs[1]
        for field_name in field_names:
            val = poly_data.GetPointData().GetArray(field_name).GetTuple1(i); 
            output += "\t" + val
    # write file
    with open(fileout, 'w') as fout:
        fout.write(header)
        fout.write(output)
    print("\tWrite ascii data file: %s" % fileout)


def InterpolateGrid(poly_data, grid_data, **kwargs):
    '''
    Inputs:
        poly_data (vtkPolyData): input data set
        grid_data (vtkGrid): grid to interpolate to
    Return:
        poly data on the new grid
    Output:
    '''
    fileout = kwargs.get('fileout', None)
    print("%s Perform interpolation onto the new grid" % Utilities.func_name())
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(poly_data)  # use the polydata
    probeFilter.SetInputData(grid_data) # Interpolate 'Source' at these points
    probeFilter.Update()
    # export to file output if a valid file is provided
    if (os.access(fileout, os.R_OK)):
        print("%s: export data" % Utilities.func_name())
        ExportPolyData(probeFilter.GetOutput(), fileout)
    return probeFilter.GetOutput()


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