# variable to substitute
# vts file:
#   PARAVIEW_FILE
# directory for images:
#   IMG_OUTPUT_DIR
# y coordinates for trench edge
#   TRENCH_EDGE_Y

class SLAB(PARAVIEW_PLOT):
    '''
    Inherit frome PARAVIEW_PLOT
    Usage is for plotting the shape of the 3D slabs
    '''
    def __init__(self, filein, **kwargs):
        PARAVIEW_PLOT.__init__(self, filein, **kwargs)
        self.eta_min = 1e19
        self.eta_max = 1e24
        

    def setup_surface_slice(self):
        '''
        Generate a visualization of a slice perpendicular to z direction, and locates
        at 10 km in depth.
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, 2000000.0, 990000.0],\
        [0.0, 0.0, 1.0], renderView=self.renderView1, name="surface_z")
        # 0: renderView.CameraPosition
        # 1: renderView.CameraFocalPoint 
        # 2: renderView.CameraParallelScale
        # 3: renderView.CameraViewUp
        _camera = [[2000000.0, 2000000.0, 14540803.753676033],\
        [2000000.0, 2000000.0, 990000.0], 4019667.7972010756, [0.0, 1.0, 0.0]]
        _color_bar = [[0.3417165169919925, 0.12791929382093295], 0.3299999999999999]
        self.camera_dict['slice_surface_z'] = _camera
        self.colorbar_dict['slice_surface_z'] = _color_bar
        adjust_camera(self.renderView1, _camera[0],_camera[1], _camera[2], _camera[3])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view
    
    def setup_trench_slice_center(self):
        '''
        Generate a visualization of a slice perpendicular to y direction, and locates at trench center (y=0)
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, 1.0, 500000.0],\
        [0.0, 1.0, 0.0], renderView=self.renderView1, name="trench_center_y")
        _camera = [[2000000.0, -14540803.753676033, 400000.0],\
        [TRENCH_INITIAL, 1.0, 400000.0],  746067.2208702065, [0.0, 0.0, 1.0]]
        _color_bar = [[0.3376097408112943, 0.03], 0.36285420944558516]
        self.camera_dict['slice_trench_center_y'] = _camera
        self.colorbar_dict['slice_trench_center_y'] = _color_bar
        adjust_camera(self.renderView1, _camera[0],_camera[1], _camera[2], _camera[3])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view

    def setup_trench_slice_edge(self):
        '''
        Generate a visualization of a slice perpendicular to y direction, and locates at trench center (y=0)
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, TRENCH_EDGE_Y, 500000.0],\
        [0.0, 1.0, 0.0], renderView=self.renderView1, name="trench_edge_y")
        _camera = [[2000000.0, -14540803.753676033, 400000.0],\
        [TRENCH_INITIAL, TRENCH_EDGE_Y, 400000.0],  746067.2208702065, [0.0, 0.0, 1.0]]
        _color_bar = [[0.3376097408112943, 0.03], 0.36285420944558516]
        self.camera_dict['slice_trench_edge_y'] = _camera
        self.colorbar_dict['slice_trench_edge_y'] = _color_bar
        adjust_camera(self.renderView1, _camera[0],_camera[1], _camera[2], _camera[3])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view
    
    def setup_slab_iso_volume_upper(self):
        '''
        Generate a visualization of the iso volume of the upper crust
        '''
        isoVolume1, isoVolume1Display, _ = add_isovolume(self.solutionpvd, "sp_lower", (0.8, 1.0), name="slab_lower")
        _camera = [[9627261.624429794, 4705025.599150724, 1231385.7667014154],\
        [TRENCH_INITIAL, 500915.187499999, 500000.00000000023], 2e6, [0.0, 0.0, 1.0]]
        self.camera_dict['isoVolume_slab_lower'] = _camera
        Hide(isoVolume1, self.renderView1)  # hide data in view

    def setup_active_clip(self):
        '''
        clip the model domain along both the x and y axis.
        '''
        # add clip 0
        clip_active0, active0Display, _ = add_clip(self.solutionpvd, [TRENCH_INITIAL + 3e6, 2e6, 5e5],\
        [1.0, 0.0, 0.0], name="active_0")
        Hide(clip_active0, self.renderView1)  # hide data in view
        # add clip 1 
        clip_active1, active1Display, _ = add_clip(clip_active0, [1e6, TRENCH_EDGE_Y + 5e5 , 5e5],\
        [0.0, 1.0, 0.0], name="active_1")
        Hide(clip_active1, self.renderView1)  # hide data in view
        # add glyph
        add_glyph("clip_active_1", "velocity", 2e5, 500)

    
    def setup_stream_tracer(self, _source):
        '''
        setup stream tracer
        '''
        # show the slice
        active_clip = FindSource(_source)
        SetActiveSource(active_clip)
        renderView1 = GetActiveViewOrCreate('RenderView')

        # create a new 'Stream Tracer'
        streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=active_clip, SeedType='Line')
        streamTracer1.Vectors = ['POINTS', 'velocity']
        streamTracer1.MaximumStreamlineLength = 7383000.0

        # init the 'Line' selected for 'SeedType'
        streamTracer1.SeedType.Point2 = [0.0, 0.0, 0.0]
        streamTracer1.SeedType.Point2 = [7383000.0, 0.0, 2890000.0]

        # show data in view
        Show(streamTracer1, renderView1, 'GeometryRepresentation')
        # streamTracer1Display = GetDisplayProperties(streamTracer1, view=renderView1)
        Hide3DWidgets()
        Hide(streamTracer1, renderView1)

        # add the 2nd one
        streamTracer2 = StreamTracer(registrationName='StreamTracer2', Input=active_clip, SeedType='Line')
        streamTracer2.Vectors = ['POINTS', 'velocity']
        streamTracer2.MaximumStreamlineLength = 7383000.0


        # init the 'sphere' selected for 'SeedType'
        streamTracer2.SeedType = 'Point Cloud'
        streamTracer2.SeedType.Radius = 300e3
        streamTracer2.SeedType.Center = [4000e3, 300e3, 2500e3]
        streamTracer2.SeedType.NumberOfPoints = 100

        Show(streamTracer2, renderView1, 'GeometryRepresentation')
        Hide3DWidgets()
        Hide(streamTracer2, renderView1)

    def plot_slice(self, _source, field_name, **kwargs):
        '''
        Plot surface slice
        After the the plots are setup, this function is then called to execute the exportation of figures.
        '''
        fig_resolution = (974, 793)
        hide_plot = kwargs.get("hide_plot", True)
        use_log = kwargs.get("use_log", False)
        _color = kwargs.get("color", None)
        lim = kwargs.get("lim", None)
        show_axis = kwargs.get("show_axis", True)
        invert_color = kwargs.get("invert_color", False)
        # show the slice
        slice1 = FindSource(_source)
        SetActiveSource(slice1)
        renderView1 = GetActiveViewOrCreate('RenderView')
        if show_axis:
            renderView1.AxesGrid.Visibility = 1 
        _display = Show(slice1, renderView1, 'GeometryRepresentation')
        slice_Display = GetDisplayProperties(slice1, view=renderView1)
        # adjust the field to plot
        ColorBy(slice_Display, ('POINTS', field_name))
        fieldLUT = GetColorTransferFunction(field_name)  # get color transfer map
        fieldPWF = GetOpacityTransferFunction(field_name) # get opacity transfer function/opacity map
        if lim != None:
            assert(len(lim)==2 and lim[0] < lim[1])
            fieldLUT.RescaleTransferFunction(lim[0], lim[1]) # Rescale transfer function
            fieldPWF.RescaleTransferFunction(lim[0], lim[1]) # Rescale transfer function
        if use_log:
            fieldLUT.UseLogScale = 1
        if _color != None:
            fieldLUT.ApplyPreset(_color, True)
        if invert_color:
            fieldLUT.InvertTransferFunction()
        # adjust colorbar and camera
        _camera = self.camera_dict[_source]
        adjust_slice_colorbar_camera(self.renderView1, fieldLUT, _camera)
        # adjust colorbar, if there is a configuration for the source
        _display.SetScalarBarVisibility(self.renderView1, True)
        try:
            _color_bar = self.colorbar_dict[_source]
            adjust_color_bar(fieldLUT, self.renderView1, _color_bar) 
        except KeyError:
            pass
        # get layout
        layout1 = GetLayout()
        # layout/tab size in pixels
        layout1.SetSize(fig_resolution[0], fig_resolution[1])
        # save figure
        file_out = os.path.join(self.output_dir, "%s_%s_%.4e.png" % (_source, field_name, self.time))
        SaveScreenshot(file_out, self.renderView1, ImageResolution=fig_resolution)
        print("Figure saved: %s" % file_out)
        if hide_plot:
            Hide(slice1, renderView1)  # hide data in view
        
    def plot_slab_volume(self, field_name, **kwargs):
        fig_resolution = (974, 793)
        lim = kwargs.get('lim', None)
        use_log = kwargs.get('use_log', False)
        _color = kwargs.get('color', None)
        invert_color = kwargs.get('invert_color', False)
        lim_slice = kwargs.get('lim_slice', None)
        use_log_slice = kwargs.get('use_log_slice', False)
        color_slice = kwargs.get('color_slice', None)
        invert_color_slice = kwargs.get('invert_color_slice', False)
        _source = "isoVolume_slab_lower"
        _source1 = "slice_trench_center_y"
        _source2 = "glyph_0"
        show_axis = True
        hide_plot = True
        field_name_slice = "viscosity"
        # show the volume
        isoVolume1 = FindSource(_source)
        SetActiveSource(isoVolume1)
        renderView1 = GetActiveViewOrCreate('RenderView')
        if show_axis:
            renderView1.AxesGrid.Visibility = 1 
        _display = Show(isoVolume1, renderView1, 'GeometryRepresentation')
        isoVolume1Display = GetDisplayProperties(isoVolume1, view=renderView1)
        # adjust the field to plot for the volume (moho surface)
        ColorBy(isoVolume1Display, ('POINTS', field_name))
        fieldLUT = GetColorTransferFunction(field_name)  # get color transfer map
        fieldPWF = GetOpacityTransferFunction(field_name) # get opacity transfer function/opacity map
        if lim != None:
            assert(len(lim)==2 and lim[0] < lim[1])
            fieldLUT.RescaleTransferFunction(lim[0], lim[1]) # Rescale transfer function
            fieldPWF.RescaleTransferFunction(lim[0], lim[1]) # Rescale transfer function
        if use_log:
            fieldLUT.UseLogScale = 1
        if _color != None:
            fieldLUT.ApplyPreset(_color, True)
        if invert_color:
            fieldLUT.InvertTransferFunction()
        # adjust colorbar for the volume, if there is a configuration for the source
        _display.SetScalarBarVisibility(renderView1, True)
        color_bar_volume = [[0.53, 0.03], 0.36285420944558516]
        adjust_color_bar(fieldLUT, renderView1, color_bar_volume) 
        # show the slice
        slice1 = FindSource(_source1)
        SetActiveSource(slice1)
        renderView1 = GetActiveViewOrCreate('RenderView')
        _display = Show(slice1, renderView1, 'GeometryRepresentation')
        slice_Display = GetDisplayProperties(slice1, view=renderView1)
        # adjust the field to plot for the slice
        ColorBy(slice_Display, ('POINTS', field_name_slice))
        fieldLUTslice = GetColorTransferFunction(field_name_slice)  # get color transfer map
        fieldPWFslice = GetOpacityTransferFunction(field_name_slice) # get opacity transfer function/opacity map
        if lim_slice != None:
            assert(len(lim_slice)==2 and lim_slice[0] < lim_slice[1])
            fieldLUTslice.RescaleTransferFunction(lim_slice[0], lim_slice[1]) # Rescale transfer function
            fieldPWFslice.RescaleTransferFunction(lim_slice[0], lim_slice[1]) # Rescale transfer function
        if use_log_slice:
            fieldLUTslice.UseLogScale = 1
        if color_slice != None:
            fieldLUTslice.ApplyPreset(color_slice, True)
        if invert_color_slice:
            fieldLUTslice.InvertTransferFunction()
        # adjust colorbar for the slice
        _display.SetScalarBarVisibility(renderView1, True)
        color_bar_slice = [[0.03, 0.03], 0.36285420944558516]  # this doesn't work for now
        adjust_color_bar(fieldLUTslice, renderView1, color_bar_slice)
        # show the glyph
        glyph1 = FindSource("Glyph1")
        SetActiveSource(glyph1)
        renderView1 = GetActiveViewOrCreate('RenderView')
        if show_axis:
            renderView1.AxesGrid.Visibility = 1 
        glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
        # adjust the colormap for the glyph
        ColorBy(glyph1Display, ('POINTS', 'velocity', 'Magnitude'))
        custom_range = [0.01, 0.1]
        velocityLUT = GetColorTransferFunction('velocity', glyph1Display, separate=True)
        velocityLUT.ApplyPreset('roma', True) # this doesn't work
        velocityLUT.InvertTransferFunction()
        if custom_range != None:
            assert(len(custom_range) == 2)
            # Rescale transfer function
            velocityLUT.RescaleTransferFunction(custom_range[0], custom_range[1])
        # adjust colorbar for the glyph
        glyph1Display.SetScalarBarVisibility(renderView1, True)
        color_bar_slice = [[0.03, 0.53], 0.36285420944558516]  # this doesn't work for now
        adjust_color_bar(velocityLUT, renderView1, color_bar_slice)
        # adjust colorbar and camera
        _camera = self.camera_dict[_source]
        adjust_camera(renderView1, _camera[0],_camera[1], _camera[2], _camera[3])
        # save figure
        file_out = os.path.join(self.output_dir, "sp_lower_%.4e.png" % (self.time))
        SaveScreenshot(file_out, renderView1, ImageResolution=fig_resolution)
        print("Figure saved: %s" % file_out)
        if hide_plot:
            Hide(isoVolume1, renderView1)  # hide data in view

    def plot_slice(self, _source): 
        '''
        plot a slice
        '''
        # get active view and source
        renderView1 = GetActiveViewOrCreate('RenderView')
        
        field1 = "T"
        field2 = "viscosity"

        # part 1: plot the center slice 
        # part 1a: temperature
        source1 = FindSource(_source)
        
        source1Display = Show(source1, renderView1, 'GeometryRepresentation')
        
        # get the original plot field, in order to hide
        # redundant color in later codes
        field0 = source1Display.ColorArrayName[1]
        field0LUT = GetColorTransferFunction(field0)

        # get color transfer function/color map for 'field'
        field1LUT = GetColorTransferFunction(field1)
        field1PWF = GetOpacityTransferFunction(field1)

        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field1, 'Magnitude'))
        HideScalarBarIfNotNeeded(field0LUT, renderView1)
        # high redundant colorbar
        source1Display.SetScalarBarVisibility(renderView1, True)
        
        # reset limit
        field1LUT.RescaleTransferFunction(273.0, 2273.0) # Rescale transfer function
        field1PWF.RescaleTransferFunction(273.0, 2273.0) # Rescale transfer function
        
        # colorbar position
        field1LUTColorBar = GetScalarBar(field1LUT, renderView1)
        field1LUTColorBar.Orientation = 'Horizontal'
        field1LUTColorBar.WindowLocation = 'Any Location'
        field1LUTColorBar.Position = [0.041, 0.908]
        field1LUTColorBar.ScalarBarLength = 0.33
        
        # Properties modified on renderView1.AxesGrid
        renderView1.AxesGrid.Visibility = 1
        # hide the orientation axes
        renderView1.OrientationAxesVisibility = 0

        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(1350, 704)
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [4688116.467077909, -15316645.446003344, 2608614.9857629854]
        renderView1.CameraFocalPoint = [4688116.467077909, 1.0, 2608614.9857629854]
        renderView1.CameraViewUp = [0.0, 0.0, 1.0]
        renderView1.CameraParallelScale = 486990.898610654

        # save figure
        fig_path = os.path.join(self.pv_output_dir, "%s_slice_center_t%.4e.pdf" % (field1, self.time))
        ExportView(fig_path, view=renderView1)

        # part 1b: viscosity
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
        field2PWF = GetOpacityTransferFunction(field2)
        
        # convert to log space
        field2LUT.MapControlPointsToLogSpace()
        field2LUT.UseLogScale = 1
        field2LUT.ApplyPreset('roma', True)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field2, 'Magnitude'))
        source1Display.SetScalarBarVisibility(renderView1, True)
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        # Hide the scalar bar for the first field color map
        HideScalarBarIfNotNeeded(field1LUT, renderView1)

        # colorbar position
        field2LUTColorBar = GetScalarBar(field2LUT, renderView1)
        field2LUTColorBar.Orientation = 'Horizontal'
        field2LUTColorBar.WindowLocation = 'Any Location'
        field2LUTColorBar.Position = [0.041, 0.908]
        field2LUTColorBar.ScalarBarLength = 0.33
        
        # Properties modified on renderView1.AxesGrid
        renderView1.AxesGrid.Visibility = 1
        # hide the orientation axes
        renderView1.OrientationAxesVisibility = 0

        # save figure
        fig_path = os.path.join(self.pv_output_dir, "%s_%s_t%.4e.pdf" % (_source, field2, self.time))
        ExportView(fig_path, view=renderView1)

        # hide plots 
        Hide(source1, renderView1)
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
    
    def plot_step(self): 
        '''
        plot a step
        '''
        self.plot_slice("slice_trench_center_y")
        self.plot_slice("slice_trench_edge_y")
        pass


def adjust_slice_colorbar_camera(renderView, colorLUT, _camera):
    '''
    adjust colorbar and camera for a slice
    Inputs:
        renderView: an instance of the rendered view
        colorLUT: an instance of the colorbar
    '''
    assert(len(_camera) == 4)
    # adjust colorbar
    colorLUTColorBar = GetScalarBar(colorLUT, renderView)
    colorLUTColorBar.WindowLocation = 'Any Location'
    colorLUTColorBar.ScalarBarLength = 0.33000000000000007
    colorLUTColorBar.Orientation = 'Horizontal'
    colorLUTColorBar.Position = [0.3458232931726907, 0.2540226986128623]
    # adjust camera
    adjust_camera(renderView, _camera[0],_camera[1], _camera[2], _camera[3])


def main():
    all_available_graphical_snapshots = ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS
    all_available_graphical_times = ALL_AVAILABLE_GRAPHICAL_TIMES
    assert(len(all_available_graphical_snapshots) == len(all_available_graphical_times))
    # First, make directory for images if it's not there
    if not os.path.isdir("IMG_OUTPUT_DIR"):
        os.mkdir("IMG_OUTPUT_DIR")
    # Process and generate plots
    Slab = SLAB("PARAVIEW_FILE", output_dir="IMG_OUTPUT_DIR")
    Slab.setup_surface_slice()
    Slab.setup_trench_slice_center()
    Slab.setup_trench_slice_edge()
    Slab.setup_slab_iso_volume_upper()
    Slab.setup_active_clip()
    Slab.setup_stream_tracer('clip_active_1')
    # First number is the number of initial adaptive refinements
    # Second one is the snapshot to plot
    # here we prefer to use a series of snapshots.
    # If this doesn't work, we will use a single snapshot
    steps = GRAPHICAL_STEPS
    if not steps == []:
        for step in steps:
            # check that snapshot is valid
            snapshot = INITIAL_ADAPTIVE_REFINEMENT+step
            if snapshot in all_available_graphical_snapshots:
                idx = all_available_graphical_snapshots.index(snapshot)
                _time =  all_available_graphical_times[idx]
                Slab.goto_time(_time)
                Slab.plot_step()
               #  Slab(INITIAL_ADAPTIVE_REFINEMENT+step)
            else:
                print ("step %s is not valid. There is no output" % step)
    else:
        snapshot = SINGLE_SNAPSHOT
        idx = all_available_graphical_snapshots.index(snapshot)
        _time =  all_available_graphical_times[idx]
        Slab.goto_time(_time)
        # Slab.plot_slice()
        # Slab(INITIAL_ADAPTIVE_REFINEMENT+SINGLE_SNAPSHOT)


main()
