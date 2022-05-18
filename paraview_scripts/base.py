# class for working with paraview
# This part contains a class (PARAVIEW_PLOT) of basic operations
# e.g. open file vtk files, add plots, choose color schemes, etc
# Usage of this class is to be combined with classes defined for
# individual visualizations

import os
# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

class PARAVIEW_PLOT():

    def __init__(self, filein, **kwargs):
        """
        initiate
        Args:
            filein(str): path of vtu file
            plots(list): lists of plot
        """
        # get variables
        self.output_dir = kwargs.get('output_dir', ".")
        # all the variables we want to plot, one in these:
        # 'velocity', 'p', 'T', 'sp_upper', 'sp_lower', 'density', 'viscosity', 
        # 'current_cohesions', 'current_friction_angles', 'plastic_yielding', 'dislocation_viscosity', 
        # 'diffusion_viscosity', 'peierls_viscosity', 'error_indicator']
        self.view_solution_pvd = False  # if the solutionpvd is included as a view
        self.all_variables = ['velocity', 'p', 'T', 'sp_upper', 'sp_lower', 'density', 'viscosity']
        # open data base
        self.filein = filein
        # assign initial values 
        self.idxs = {}
        self.all_idxs = 0
        self.time_snap = 0
        self.solutionpvd = PVDReader(registrationName='solution.pvd', FileName=filein)
        self.solutionpvd.PointArrays = self.all_variables
        # get animation scene
        animationScene1 = GetAnimationScene()
        # update animation scene based on data timesteps
        animationScene1.UpdateAnimationUsingDataTimeSteps()
        # get active view
        self.renderView1 = GetActiveViewOrCreate('RenderView')
        if self.view_solution_pvd:
            # show data in view
            solutionpvdDisplay = Show(self.solutionpvd, self.renderView1, 'UnstructuredGridRepresentation')
            # hide data in view
            Hide(self.solutionpvd, self.renderView1)
        # update the view to ensure updated data information
        self.renderView1.Update()




def add_slice(solutionpvd, field, Origin, Normal, **kwargs):
    '''
    create a new 'Slice'
    Inputs:
        solutionpvd: a solutionpvd object of paraview
        field: field to show
        Origin: the origin point of visualization
        Normal: the normal direction of the slice
    '''
    # get additional variables
    renderView0 = kwargs.get('renderView', None)
    # create slice
    slice1 = Slice(registrationName='Slice1', Input=solutionpvd)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = Origin
    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [2000000.0, 2000000.0, 500000.0]
    # Modify the normal direction
    slice1.SliceType.Normal = Normal
    # get active view
    if renderView0 == None:
        renderView1 = GetActiveViewOrCreate('RenderView')
    else:
        renderView1 = renderView0
    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    adjust_slice_color(slice1Display, field, renderView1)
    renderView1.Update()
    return slice1, slice1Display, renderView1

    
def adjust_slice_color(slice1Display, field, renderView):
    '''
    adjust the color scheme of a slice filter
    Inputs:
        slice1Display : a display object we created before this
        field : field to show
        renderView: a renderView object from paraview
    '''
    # set scalar coloring
    ColorBy(slice1Display, ('POINTS', field))
    # Hide the scalar bar for this color map if no visible data is colored by it.
    # HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)
    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView, True)
    # get color transfer function/color map for field
    sp_upperLUT = GetColorTransferFunction(field)
    # get opacity transfer function/opacity map for field
    sp_upperPWF = GetOpacityTransferFunction(field)


def adjust_camera(renderView, position, focalPoint, parallelScale, viewUp):
    '''
    adjust camera
    Inputs:
        renderView: a "renderView" object of Paraview
    '''
    # reset view to fit data
    renderView.ResetCamera(False)
    # current camera placement for renderView
    renderView.InteractionMode = '2D'
    renderView.CameraPosition = position
    renderView.CameraFocalPoint = focalPoint
    renderView.CameraParallelScale = parallelScale
    renderView.CameraViewUp = viewUp