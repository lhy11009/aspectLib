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
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
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
        # For splitting the domain
        self.spacing_n = 0
        self.spacing = None
        self.spacing_cell_ids = None
        # for storing an array of cells
        self.cell_array = None
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

    def ReadFile(self, filein, **kwargs):
        '''
        Read file
        Inputs:
            filein (str): path to a input file
        '''
        quiet = kwargs.get("quiet", False)
        if not quiet:
            print("%s started" % Utilities.func_name())
        start = time.time()
        assert(os.path.isfile(filein))
        file_extension = filein.split('.')[-1]
        if file_extension == 'pvtu':
            self.reader = vtk.vtkXMLPUnstructuredGridReader()
        elif file_extension == 'vtu':
            self.reader = vtk.vtkXMLUnstructuredGridReader()
        else:
            raise TypeError("%s: Wrong type of file" % Utilities.func_name())
        self.reader.SetFileName(filein)
        self.reader.Update()
        end = time.time()
        if not quiet:
            print("\tReadFile: %s" % filein)
            print("\t%s, takes %.2f s" % (Utilities.func_name(), end - start))
    
    def ConstructPolyData(self, field_names, **kwargs):
        '''
        construct poly data
        Inputs:
            field_names (list of field names)
            kwargs:
                include_cell_center - include_cell_center in the poly_data
        '''
        quiet = kwargs.get("quiet", False)
        if not quiet:
            print("%s started" % Utilities.func_name())
        start = time.time()
        include_cell_center = kwargs.get('include_cell_center', False)
        fix_cell_value = kwargs.get('fix_cell_value', True)  # fix the value of cell centers
        construct_Tdiff = kwargs.get("construct_Tdiff", False)
        assert(type(field_names) == list and len(field_names) > 0)
        noP = 0
        noC = 0
        grid = self.reader.GetOutput()
        data_set = self.reader.GetOutputAsDataSet()
        points = grid.GetPoints()
        cells = grid.GetCells()
        self.cells = cells
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
        # include cell centers
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
        message = "\tConstructPolyData: %d * (%d + %d) entries in the polydata imported and %d * (%d + %d) points in the data at cell center. \
import data takes %f, interpolate cell center data takes %f"\
        % (noP, self.dim, len(field_names), noC, self.dim, len(field_names), time_import - start, time_center - time_import)
        if not quiet:
            print(message)
        end = time.time()
        if not quiet:
            print("\tConstruct polydata, takes %.2f s" % (end - start))

    def SplitInSpace(self, spacing, **kwargs):
        '''
        Split the space
        '''
        print("%s started" % Utilities.func_name())
        start = time.time()
        # options
        geometry = kwargs.get('geometry', 'box')
        if geometry == "box":
            is_cartesian = True
        elif geometry == "chunk":
            is_cartesian = False
        else:
            raise ValueError("%s: geometry needs to be cartesian or chunk" % Utilities.func_name())
        dim = kwargs.get('dim', 2)

        # Get the bounds
        domain_bounds = self.i_poly_data.GetPoints().GetBounds()

        # spacing: split the domain into smaller space, record the number of the 3 axis and the total number
        self.spacing = spacing
        if dim == 2:
            return NotImplementedError
        elif dim == 3:
            spacing_x, spacing_y, spacing_z = spacing[0], spacing[1], spacing[2]
            interval_x = (domain_bounds[1] - domain_bounds[0]) / spacing_x
            interval_y = (domain_bounds[3] - domain_bounds[2]) / spacing_y
            interval_z = (domain_bounds[5] - domain_bounds[4]) / spacing_z
            spacing_n = spacing_x * spacing_y * spacing_z
            for ix in range(spacing_x):
                for jy in range(spacing_y):
                    for kz in range(spacing_z):
                        spacing_idx = kz + jy * spacing_z + ix * spacing_y * spacing_z
        self.spacing_n = spacing_n

        # distribute cells by looking up the cell center
        self.spacing_cell_ids = [[] for i in range(self.spacing_n)]
        centers = vtk.vtkCellCenters()  # test the utilities of centers
        centers.SetInputData(self.i_poly_data)
        centers.Update()
        cell_centers = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
        for iC in range(cell_centers.shape[0]):
            cell_coordinates = cell_centers[iC]
            cell_center_x, cell_center_y, cell_center_z = cell_coordinates[0], cell_coordinates[1], cell_coordinates[2]
            ix = int((cell_center_x - domain_bounds[0]) // interval_x)
            jy = int((cell_center_y - domain_bounds[2]) // interval_y)
            kz = int((cell_center_z - domain_bounds[4]) // interval_z)
            spacing_idx = kz + jy * spacing_z + ix * spacing_y * spacing_z
            self.spacing_cell_ids[spacing_idx].append(iC)
        end = time.time()
        print("\t%s takes %.2f s" % (Utilities.func_name(), end-start))
    
    def InterpolateDomain(self, points, fields, **kwargs):
        '''
        Run interpolation
        '''
        quiet = kwargs.get("quiet", True)
        if not quiet:
            print("%s started" % Utilities.func_name())
        start = time.time()
        cells_vtk = kwargs.get("cells_vtk", None) # for set connectivity
        points_found = kwargs.get("points_found", None)
        output_poly_data = kwargs.get("output_poly_data", True)
        interpolated_data = kwargs.get("interpolated_data", None)
        apply_additional_chunk_search = kwargs.get('apply_additional_chunk_search', True)
        is_box = (self.geometry == "box")
        if not is_box:
            Utilities.my_assert(self.geometry == "chunk", ValueError, "%s: we only handle box and chunk geometry" % Utilities.func_name())
        # check point dimension
        if points.shape[1] == 2:
            pass
        elif points.shape[1] == 3:
            pass
        else:
            raise ValueError("%s: points.shape[1] needs to either 2 or 3" % Utilities.func_name())
        # initiate a ndarray to store found information 
        if points_found is None:
            points_found = np.zeros(points.shape[0], dtype=int)
        # Get the bounds
        domain_bounds = self.i_poly_data.GetPoints().GetBounds()
        # prepare the datafield to interpolate
        # field_type: 0 - scalar, 1 - vector
        raw_data = []
        field_type = []
        n_vector = 0
        for field in fields:
            if field in ["velocity"]:
                raw_data.append(self.i_poly_data.GetPointData().GetVectors(field))
                field_type.append(1)
                n_vector += 1
            else:
                raw_data.append(self.i_poly_data.GetPointData().GetArray(field))
                field_type.append(0)
        if interpolated_data is None:
            interpolated_data = np.zeros((len(fields), points.shape[0]))
        if n_vector > 0:
            interpolated_vector = [ [[0.0, 0.0, 0.0] for j in range(points.shape[0])] for i in range(n_vector)]
        end = time.time()
        if not quiet:
            print("\tInitiating takes %2f s" % (end - start))
            if (not is_box) and apply_additional_chunk_search:
                print("\tApply additional chunk search")
        
        # loop over the points, find cells containing the points and interpolate from the cell points
        start = end
        points_in_cell = [[0.0, 0.0, 0.0] for i in range(4)] # for 2d, parse vtk cell, chunk case
        n_found = 0
        n_out_of_bound = 0
        n_not_found = 0
        for i in range(points.shape[0]):
            if points_found[i]:
                # skip points found in other files
                continue
            if not PointInBound2D(points[i], domain_bounds):
                # coordinates out of range, append invalid values
                n_out_of_bound += 1
                continue
            for iC in range(self.i_poly_data.GetNumberOfCells()):
                cell = self.i_poly_data.GetCell(iC)
                bound = cell.GetBounds()
                if PointInBound2D(points[i], bound):
                    # mark points inside the cell iC and mark found of point i
                    # box: simply use the cell bound
                    # chunk: check first the cell bound, then convert both the query point and the bound
                    # to spherical coordinates and check again
                    if is_box:
                        found = True
                        points_found[i] = 1
                        cell_found = cell
                        break
                    else:
                        if apply_additional_chunk_search:
                            r, theta, phi = Utilities.cart2sph(points[i][0], points[i][1], points[i][2])
                            cell_points = cell.GetPoints()
                            for i_pc in range(cell.GetNumberOfPoints()):
                                point_id = cell.GetPointId(i_pc)
                                cell_points.GetPoint(i_pc, points_in_cell[i_pc])
                            sph_bounds_cell = Utilities.SphBound(points_in_cell)
                            if PointInBound([phi, theta, r], sph_bounds_cell):
                                found = True
                                points_found[i] = 1
                                cell_found = cell
                                break
                        else:
                            found = True
                            points_found[i] = 1
                            cell_found = cell
                            break
            if found:
                n_found += 1
                # Prepare variables for EvaluatePosition
                closest_point = [0.0, 0.0, 0.0]
                sub_id = vtk.reference(0)
                dist2 = vtk.reference(0.0)
                pcoords = [0.0, 0.0, 0.0]
                weights = [0.0] * cell.GetNumberOfPoints()
                # Evaluate the position to check if the point is inside the cell
                inside = cell.EvaluatePosition(points[i], closest_point, sub_id, pcoords, dist2, weights)
                for i_f in range(len(fields)):
                    fdata = raw_data[i_f]
                    if field_type[i_f] == 0:
                        interpolated_val = 0.0
                        for i_pc in range(cell_found.GetNumberOfPoints()):
                            point_id = cell_found.GetPointId(i_pc)
                            value = fdata.GetTuple1(point_id)  # Assuming scalar data
                            interpolated_val += value * weights[i_pc]
                        interpolated_data[i_f][i] = interpolated_val
                    elif field_type[i_f] == 1:
                        interpolated_val = np.array([0.0, 0.0, 0.0])
                        for i_pc in range(cell_found.GetNumberOfPoints()):
                            point_id = cell_found.GetPointId(i_pc)
                            value = fdata.GetTuple(point_id)  # Assuming scalar data
                            interpolated_val += np.array(value) * weights[i_pc]
                        interpolated_vector[0][i] = interpolated_val
            else:
                n_not_found += 1 
        total_n_found = np.sum(points_found==1)
        end = time.time()
        if not quiet:
            print("\t%s Searched %d points (%d found total; %d found current file; %d out of bound; %d not found)" % (Utilities.func_name(), points.shape[0], total_n_found, n_found, n_out_of_bound, n_not_found))
            print("\tSearching takes %2f s" % (end - start))

        # construct new polydata
        # We construct the polydata ourself, now it only works for field data
        if output_poly_data:
            o_poly_data = vtk.vtkPolyData()
            points_vtk = vtk.vtkPoints()
            for i in range(points.shape[0]):
                points_vtk.InsertNextPoint(points[i])
            o_poly_data.SetPoints(points_vtk) # insert points
            if cells_vtk is not None:
                o_poly_data.SetPolys(cells_vtk)
            for i_f in range(len(fields)):
                interpolated_array = numpy_to_vtk(interpolated_data[i_f])
                interpolated_array.SetName(fields[i_f])
                if i_f == 0:
                    # the first array
                    o_poly_data.GetPointData().SetScalars(interpolated_array)
                else:
                    # following arrays
                    o_poly_data.GetPointData().AddArray(interpolated_array)
        else:
            o_poly_data = None

        return o_poly_data, points_found, interpolated_data, interpolated_vector

    
    def InterpolateSplitSpacing(self, points, fields, **kwargs):
        '''
        Run interpolation from split spacing
        '''
        print("%s started" % Utilities.func_name())
        start = time.time()
        assert(self.spacing_n > 0)
        cells_vtk = kwargs.get("cells_vtk", None) # for set connectivity
        split_perturbation = kwargs.get("split_perturbation", 1) # number of split to acess during interpolation
        debug = kwargs.get("debug", False) # print debug information
        points_found = kwargs.get("points_found", None)
        output_poly_data = kwargs.get("output_poly_data", True)
        interpolated_data = kwargs.get("interpolated_data", None)
        apply_additional_chunk_search = kwargs.get('apply_additional_chunk_search', True)
        is_box = (self.geometry == "box")
        if not is_box:
            Utilities.my_assert(self.geometry == "chunk", ValueError, "%s: we only handle box and chunk geometry" % Utilities.func_name())
        # check point dimension
        if points.shape[1] == 2:
            raise NotImplementedError()
        elif points.shape[1] == 3:
            pass
        else:
            raise ValueError("%s: points.shape[1] needs to either 2 or 3" % Utilities.func_name())
        # initiate a ndarray to store found information 
        if points_found is None:
            points_found = np.zeros(points.shape[0], dtype=int)
        # Get the bounds
        domain_bounds = self.i_poly_data.GetPoints().GetBounds()
        # Extract the range for x, y, and z
        spacing_x, spacing_y, spacing_z = self.spacing[0], self.spacing[1], self.spacing[2]
        interval_x = (domain_bounds[1] - domain_bounds[0]) / spacing_x
        interval_y = (domain_bounds[3] - domain_bounds[2]) / spacing_y
        interval_z = (domain_bounds[5] - domain_bounds[4]) / spacing_z

        # prepare the datafield to interpolate
        raw_data = []
        for field in fields:
            raw_data.append(self.i_poly_data.GetPointData().GetArray(field))
        if interpolated_data is None:
            interpolated_data = np.zeros((len(fields), points.shape[0]))

        # interpolate
        # Method to use: simply looking at the bound of each cell
        # Note there is a method of EvaluatePosition to decide wheter a point is in a cell.
        # But that doesn't work well. It give different results for the same point and cell, shame.
        # looking at the bound works well for cartesian geometry but needs double-check the spherical bound in chunk geometry
        n_found = 0  # record whether interpolation is successful
        n_out_of_bound = 0
        n_not_found = 0
        # apply search over multiple slices by perturbation on the
        # slice index in all 3 directions
        diffs = []
        for pert_x in range(split_perturbation):
            for pert_y in range(split_perturbation):
                for pert_z in range(split_perturbation):
                    diffs.append([pert_x, pert_y, pert_z])
                    if pert_x > 0:
                        diffs.append([-pert_x, pert_y, pert_z])
                        if pert_y > 0:
                            diffs.append([-pert_x, -pert_y, pert_z])
                            if pert_z > 0:
                                diffs.append([-pert_x, -pert_y, -pert_z])
                        if pert_z > 0:
                            diffs.append([-pert_x, pert_y, -pert_z])
                    if pert_y > 0:
                        diffs.append([pert_x, -pert_y, pert_z])
                        if pert_z > 0:
                            diffs.append([pert_x, -pert_y, -pert_z])
                    if pert_z > 0:
                        diffs.append([pert_x, pert_y, -pert_z])
        end = time.time()
        print("\tInitiating takes %2f s" % (end - start))
        if (not is_box) and apply_additional_chunk_search:
            print("\tApply additional chunk search")
        # loop over the points, find cells containing the points and interpolate from the cell points
        start = end
        points_in_cell = [[0.0, 0.0, 0.0] for i in range(8)] # for 3d, parse vtk cell, chunk case
        for i in range(points.shape[0]):
            if points_found[i]:
                # skip points found in other files
                continue
            if not PointInBound(points[i], domain_bounds):
                # coordinates out of range, append invalid values
                n_out_of_bound += 1
                continue
            ix = int((points[i][0] - domain_bounds[0]) // interval_x)
            jy = int((points[i][1] - domain_bounds[2]) // interval_y)
            kz = int((points[i][2] - domain_bounds[4]) // interval_z)
            # Evaluate the position
            found = False
            for diff in diffs:
                spacing_idx = Utilities.clamp(kz + diff[2], 0, spacing_z-1)\
                             + Utilities.clamp(jy + diff[1], 0, spacing_y-1) * spacing_z\
                             + Utilities.clamp(ix + diff[0], 0, spacing_x-1) * spacing_y * spacing_z
                assert(spacing_idx < self.spacing_n) # make sure we have a valid value
                for iC in self.spacing_cell_ids[spacing_idx]:
                    cell = self.i_poly_data.GetCell(iC)
                    bound = cell.GetBounds()
                    if PointInBound(points[i], bound):
                        # mark points inside the cell iC and mark found of point i
                        # box: simply use the cell bound
                        # chunk: check first the cell bound, then convert both the query point and the bound
                        # to spherical coordinates and check again
                        if is_box:
                            found = True
                            points_found[i] = 1
                            break
                        else:
                            if apply_additional_chunk_search:
                                r, theta, phi = Utilities.cart2sph(points[i][0], points[i][1], points[i][2])
                                cell_points = cell.GetPoints()
                                for i_pc in range(cell.GetNumberOfPoints()):
                                    point_id = cell.GetPointId(i_pc)
                                    cell_points.GetPoint(i_pc, points_in_cell[i_pc])
                                sph_bounds_cell = Utilities.SphBound(points_in_cell)
                                if PointInBound([phi, theta, r], sph_bounds_cell):
                                    found = True
                                    points_found[i] = 1
                                    break
                            else:
                                found = True
                                points_found[i] = 1
                                break
                if found:
                    cell_found = cell
                    break
            if found:
                n_found += 1
                # Prepare variables for EvaluatePosition
                closest_point = [0.0, 0.0, 0.0]
                sub_id = vtk.reference(0)
                dist2 = vtk.reference(0.0)
                pcoords = [0.0, 0.0, 0.0]
                weights = [0.0] * cell.GetNumberOfPoints()
                # Evaluate the position to check if the point is inside the cell
                inside = cell.EvaluatePosition(points[i], closest_point, sub_id, pcoords, dist2, weights)
                for i_f in range(len(fields)):
                    fdata = raw_data[i_f]
                    interpolated_val = 0.0
                    for i_pc in range(cell_found.GetNumberOfPoints()):
                        point_id = cell_found.GetPointId(i_pc)
                        value = fdata.GetTuple1(point_id)  # Assuming scalar data
                        interpolated_val += value * weights[i_pc]
                    interpolated_data[i_f][i] = interpolated_val
            else:
                n_not_found += 1
        total_n_found = np.sum(points_found==1)
        end = time.time()
        print("\t%s Searched %d points (%d found total; %d found current file; %d out of bound; %d not found)" % (Utilities.func_name(), points.shape[0], total_n_found, n_found, n_out_of_bound, n_not_found))
        print("\tSearching takes %2f s" % (end - start))
        # construct new polydata
        # We construct the polydata ourself, now it only works for field data
        if output_poly_data:
            o_poly_data = vtk.vtkPolyData()
            points_vtk = vtk.vtkPoints()
            for i in range(points.shape[0]):
                points_vtk.InsertNextPoint(points[i])
            o_poly_data.SetPoints(points_vtk) # insert points
            if cells_vtk is not None:
                o_poly_data.SetPolys(cells_vtk)
            for i_f in range(len(fields)):
                interpolated_array = numpy_to_vtk(interpolated_data[i_f])
                interpolated_array.SetName(fields[i_f])
                if i_f == 0:
                    # the first array
                    o_poly_data.GetPointData().SetScalars(interpolated_array)
                else:
                    # following arrays
                    o_poly_data.GetPointData().AddArray(interpolated_array)
        else:
            o_poly_data = None

        return o_poly_data, points_found, interpolated_data
            


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


def InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, **kwargs):
    '''
    Interpolation of vtu
    '''
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Ri = Visit_Options.options['INNER_RADIUS']
    Xmax = Visit_Options.options['XMAX']
    # read file
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    # construct poly data
    VtkP.ConstructPolyData(fields, include_cell_center=False)
    
    ### split the space
    VtkP.SplitInSpace(spacing, dim=3)
    
    ### Interpolate onto a new mesh
    kwargs["fields"] = fields
    o_poly_data, points_found, interpolated_data = VtkP.InterpolateSplitSpacing(target_points_np, **kwargs)
    return o_poly_data, points_found, interpolated_data

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


def PointInBound(point, bound):
    '''
    Determine whether a 3d point is within a bound
    Inputs:
        points: 3d, list or np.ndarray
        bound: 6 component list, x_min, x_max, y_min, y_max, z_min, z_min
    '''
    if type(point) == np.ndarray:
        assert(point.size == 3)
    elif type(point) == list:
        assert(len(point) == 3)
    else:
        raise TypeError
    return (point[0] >= bound[0]) and (point[0] <= bound[1]) and (point[1] >= bound[2])\
        and (point[1] <= bound[3]) and (point[2] >= bound[4]) and (point[2] <= bound[5])


def PointInBound2D(point, bound):
    '''
    Determine whether a 2d point is within a bound
    Inputs:
        points: 2d, list or np.ndarray, not we still need 3 entries to be consistent with vtk
        bound: 4 component list, x_min, x_max, y_min, y_max
    '''
    if type(point) == np.ndarray:
        assert(point.size == 3)
    elif type(point) == list:
        assert(len(point) == 3)
    else:
        raise TypeError
    return (point[0] >= bound[0]) and (point[0] <= bound[1]) and (point[1] >= bound[2])\
        and (point[1] <= bound[3])


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
    polys = kwargs.get('polys', None)
    assert(points.ndim == 2)
    if points.shape[1] == 2:
        is_2d = True
    elif points.shape[1] == 3:
        is_2d = False
    else:
        raise ValueError("%s: points.shape[1] needs to either 2 or 3" % Utilities.func_name())
    if not quiet:
        print("%s: Perform interpolation onto the new grid" % Utilities.func_name())
    start = time.time()
    grid_points = vtk.vtkPoints()
    # for point in points:
    #    grid_points.InsertNextPoint(point)
    for i in range(points.shape[0]):
        x = points[i][0]
        y = points[i][1]
        if is_2d:
            z = 0.0
        else:
            z = points[i][2]
        grid_points.InsertNextPoint(x, y, z)  # this is always x, y, z
    grid_data = vtk.vtkPolyData()
    grid_data.SetPoints(grid_points)
    end = time.time()
    if not quiet:
        print("%s: Construct vtkPoints, take %.2f s" % (Utilities.func_name(), end - start))
    start = end
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(poly_data)  # use the polydata
    probeFilter.SetInputData(grid_data) # Interpolate 'Source' at these points
    probeFilter.Update()
    o_poly_data = probeFilter.GetOutput()
    if polys is not None:
        o_poly_data.SetPolys(polys)
    end = time.time()
    if not quiet:
        print("%s: Interpolate onto new grid, take %.2f s" % (Utilities.func_name(), end - start))
    start = end
    # export to file output if a valid file is provided
    if fileout != None:
        ExportPolyData(o_poly_data, fileout, indent=4)
        end = time.time()
        print("%s: Export data, takes %.2f s" % (Utilities.func_name(), end - start))
    return o_poly_data


def ProjectVelocity(x_query, y_query, vs, geometry):
    '''
    project the velocity from vx, vy to vr, vh in different geometries
    Inputs:
        x_query - x coordinate of the point
        y_query - y coordinate of the point
        vs - (vx, vy)
        geometry - type of geometry
    ''' 
    if geometry == "chunk":
        cos_sp = x_query / (x_query**2.0 + y_query**2.0)**0.5
        sin_sp = y_query / (x_query**2.0 + y_query**2.0)**0.5
        if vs.ndim == 1:
            v_h = vs[0] * (-sin_sp) + vs[1] * cos_sp
            v_r = vs[0] * cos_sp + vs[1] * sin_sp
        elif vs.ndim == 2:
            v_h = vs[:, 0] * (-sin_sp) + vs[:, 1] * cos_sp
            v_r = vs[:, 0] * cos_sp + vs[:, 1] * sin_sp
        else:
            NotImplementedError()
    elif geometry == "box":
        if vs.ndim == 1:
            v_h = vs[0]
            v_r = vs[1]
        elif vs.ndim == 2:
            v_h = vs[:, 0]
            v_r = vs[:, 1]
        else:
            NotImplementedError()
    else:
        raise NotImplementedError()
    return v_h, v_r


def MakeTargetMesh(Visit_Options, n0, n1, d_lateral):
    '''
    Make a target mesh for slicing 3d dataset
    Inputs:
        Visit_Options - a VISIT_OPTIONS class
        n0, n1 - number of points along the 1st and 3rd dimention
        interval - this determines the interval of the slices
        d_lateral - the lateral distance, along the 2nd dimention
            take a minimum value of 1.0 to assure successful slicing of the geometry
    '''
    # get the options 
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Ri = Visit_Options.options['INNER_RADIUS']
    Xmax = Visit_Options.options['XMAX']
    N = n0 * n1
    # new mesh
    target_points_np = np.zeros((N, 3))
    if geometry == "box":
        for i0 in range(n0):
            for j1 in range(n1):
                ii = i0 * n1 + j1
                target_points_np[ii, 0] = Xmax * i0 / (n0 - 1)
                target_points_np[ii, 1] = d_lateral
                target_points_np[ii, 2] = Ro * j1 / (n1 - 1) # Ro and depth are the same in this case
    elif geometry == "chunk":
        for i0 in range(n0):
            for j1 in range(n1):
                # note we use theta = 0.0 here, but z0 = small value, this is to ensure a successful slice
                # of the chunk geometry
                # we take the slice at z = d_lateral, then x and y are between Ri and a R1 value
                if d_lateral < 1e3:
                    R1 = Ro
                else:
                    R1 = (Ro**2.0 - d_lateral**2.0)**0.5
                ii = i0 * n1 + j1
                phi = i0 / n0 * Xmax * np.pi / 180.0
                r = R1 * (j1 - 0.0) / (n1 - 0.0) + Ri * (j1 - n1)/ (0.0 - n1)  
                slice_x, slice_y, _  = Utilities.ggr2cart(0.0, phi, r) 
                target_points_np[ii, 0] = slice_x
                target_points_np[ii, 1] = slice_y
                target_points_np[ii, 2] = d_lateral
    return target_points_np


def GetVtkCells2d(n_x, n_z):
    '''
    Get a vtk cells from number of points along x, z in 2 dimensions
    '''
  
    cells = vtk.vtkCellArray()
  
    for ix in range(n_x-1):
        for jz in range(n_z-1):
          cells.InsertNextCell(4)
          # cell = vtk.vtkIdType()
          cells.InsertCellPoint(ix*n_z + jz+1)
          cells.InsertCellPoint(ix*n_z + jz)
          cells.InsertCellPoint((ix+1)*n_z + jz)
          cells.InsertCellPoint((ix+1)*n_z + jz+1)
  
    return cells


def Interpolate3dVtkCaseBeta(case_dir, VISIT_OPTIONS, vtu_snapshot, fields, mesh_options, **kwargs):
    '''
    Inputs:
        case_dir - case directory
        VISIT_OPTIONS - a matching class containing the options
        vtu_snapshot (int) - vtu snapshot
        fields (list of str) - the fields to output
        mesh_options - dict for mesh options
            type - type of meshes
            resolution - mesh resolutions
            d_lateral - if type is slice_2nd, this is the distance measured on the 2nd dimension
        kwargs - dict
            by_part - if the interpolation is performed on individua vtu file (False: on the pvtu file)
            spacing - spacing of the domain, this determines the number of slices the algorithm produce
            split_perturbation - This determines the number of slices the algorithm searches for a query point
    ####
    Note on the trade off between spacing and the split_perturbation parameters:
    # The spacing parameter tends to divide the domain into multiple spaces and accelerate
    # the process of interpolation, but make it harder to find the cell for a query point.
    # The problem is caused by the location of the cell center. When the cell is big, the cell center
    # might be far from a point within the cell, and one cell could be splited into different pieces of spacing.
    # This tackled by a large number of split perturbation, which will tell the interpolation algorithm to look into mulitple pices of spacing
    # rather than one and increases the possibility to locate the query point in a cell.
    # In application, first start with smaller spacing numbers. If the interpolation is slow, increase
    # this number.
    '''
    #options
    # case directory
    case_dir = '/mnt/lochy/ASPECT_DATA/ThDSubduction/chunk_test/chunk_initial9'
    assert(os.path.isdir(case_dir))
    # algorithm
    by_part = kwargs.get("by_part", False)
    apply_additional_chunk_search = kwargs.get("apply_additional_chunk_search", True)
    spacing = kwargs.get("spacing", [10, 10, 10])
    split_perturbation = kwargs.get("split_perturbation", 2)
    fields = ["T", "density"]
    # mesh
    _type = mesh_options['type']
    resolutions = mesh_options['resolution']
    n0 = resolutions[0]
    n1 = resolutions[1]
    if _type == "slice_2nd":
        d_lateral = mesh_options["d_lateral"]
    else:
        raise NotImplementedError

    #Initiation 
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    
    ### make a target mesh as well as a poly (connectivity in cells)
    start = time.time()
    print("Make a target mesh")
    target_points_np = MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    target_cells_vtk = GetVtkCells2d(n0, n1)
    end = time.time()
    print("\tPoints in target: %d" % (target_points_np.shape[0]))
    print("\tOperation takes %.2f s" % (end - start))
    
    ### Perform interpolation
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    print("\tTime = %.4e" % (float(_time) / 1e6))
    interpolated_data = np.zeros((len(fields), target_points_np.shape[0]))
    if by_part:
        for part in range(16):
            filein = os.path.join(case_dir, "output", "solution", "solution-%05d.%04d.vtu" % (vtu_snapshot, part))
            if part == 0:
                points_found = None
            print("-"*20 + "split" + "-"*20) # print a spliting
            print(filein)
            _, points_found, interpolated_data = InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, points_found=points_found,\
                                                        split_perturbation=split_perturbation, interpolated_data=interpolated_data, output_poly_data=False,
                                                        apply_additional_chunk_search=apply_additional_chunk_search)
            if np.sum(points_found == 1) == points_found.size:
                print("All points have been found, exiting")
                break
    else:
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
        _, points_found, interpolated_data = InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, split_perturbation=split_perturbation,\
                                                    interpolated_data=None, output_poly_data=False, apply_additional_chunk_search=apply_additional_chunk_search)
    pass

    # organize output
    o_poly_data = vtk.vtkPolyData()
    points_vtk = vtk.vtkPoints()
    for i in range(target_points_np.shape[0]):
        points_vtk.InsertNextPoint(target_points_np[i])
    o_poly_data.SetPoints(points_vtk) # insert points
    if target_cells_vtk is not None:
        o_poly_data.SetPolys(target_cells_vtk)
    for i_f in range(len(fields)):
        interpolated_array = numpy_to_vtk(interpolated_data[i_f])
        interpolated_array.SetName(fields[i_f])
        if i_f == 0:
            # the first array
            o_poly_data.GetPointData().SetScalars(interpolated_array)
        else:
            # following arrays
            o_poly_data.GetPointData().AddArray(interpolated_array)

    # write output
    dirout = os.path.join(case_dir, "vtk_outputs")
    if not os.path.isdir(dirout):
        os.mkdir(dirout)
    if _type == "slice_2nd":
        filename = "slice_2nd_%05d.vtp" % vtu_snapshot
    fileout = os.path.join(dirout, filename)
    ExportPolyData(o_poly_data, fileout, indent=4)
    assert(os.path.isfile(fileout))
    print("%s: Write output %s" % (Utilities.func_name(), fileout))
    

#------------ Utility functions ---------------- #

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