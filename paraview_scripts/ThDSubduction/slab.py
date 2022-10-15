# variable to substitute
# vts file:
#   PARAVIEW_FILE
# directory for images:
#   IMG_OUTPUT_DIR
# y coordinates for trench edge
#   TRENCH_EDGE_Y
# minimum viscosity
#   ETA_MIN
# maximum viscosity
#   ETA_MAX

class SLAB(PARAVIEW_PLOT):
    '''
    Inherit frome PARAVIEW_PLOT
    Usage is for plotting the shape of the 3D slabs
    '''
    def __init__(self, filein, **kwargs):
        PARAVIEW_PLOT.__init__(self, filein, **kwargs)
        self.eta_min = ETA_MIN
        self.eta_max = ETA_MAX
        

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
        # adjust colorbar and camera
        _camera = self.camera_dict[_source]
        adjust_camera(renderView1, _camera[0],_camera[1], _camera[2], _camera[3])
        # save figure
        file_out = os.path.join(self.output_dir, "sp_lower_%.4e.png" % (self.time))
        SaveScreenshot(file_out, renderView1, ImageResolution=fig_resolution)
        print("Figure saved: %s" % file_out)
        if hide_plot:
            Hide(isoVolume1, renderView1)  # hide data in view

    def plot_step(self): 
        '''
        plot a step
        '''
        self.plot_slice("slice_trench_center_y", "sp_upper", color='Reds', invert_color=True, lim=[0.0, 1.0])
        self.plot_slice("slice_trench_center_y", "T", color='vik')
        self.plot_slice("slice_trench_edge_y", "sp_upper", color='Reds', invert_color=True, lim=[0.0, 1.0])
        self.plot_slice("slice_surface_z", "sp_upper", color='Reds', invert_color=True, lim=[0.0, 1.0])
        self.plot_slice("slice_trench_center_y", "viscosity", use_log=True, color='roma', lim=[self.eta_min, self.eta_max])
        self.plot_slab_volume("T", color='vik', lim_slice=[self.eta_min, self.eta_max], use_log_slice=True, color_slice='roma')


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
