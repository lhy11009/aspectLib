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
# rotation of the domain
#   ROTATION_ANGLE
# geometry of the model
#   GEOMETRY
# plot model axis
#   PLOT_AXIS
# types of plot (a list, e.g. upper mante / whole)
#   PLOT_TYPES

class SLAB(PARAVIEW_PLOT):
    '''
    Inherit frome PARAVIEW_PLOT
    Usage is for plotting the shape of the 3D slabs
    '''
    def __init__(self, filein, **kwargs):
        '''
        Initiation
        '''
        kwargs['project'] = "TwoDSubduction"
        PARAVIEW_PLOT.__init__(self, filein, **kwargs)
        self.eta_min = ETA_MIN
        self.eta_max = ETA_MAX
        self.T_min = 273.0
        self.T_max = 2273.0
        self.camera_dict['twod_upper_mantle'] = [[0.0, 6e6, 2.5e7],[0.0, 6e6, 0.0], 5.4e5, None]
        apply_rotation("solution.pvd", [0.0, 0.0, 0.0], [0.0, 0.0, ROTATION_ANGLE], registrationName="Transform1")
        add_plot("Transform1", "viscosity", use_log=True, lim=[self.eta_min, self.eta_max], color="roma")
        add_plot("Transform1", "T", lim=[self.T_min, self.T_max], color="vik")
        add_glyph1("Transform1", "velocity", 1e6, registrationName="Glyph1")
        add_deformation_mechanism("Transform1", registrationName="pFilter_DM")


    def plot_step(self, **kwargs): 
        '''
        plot a step
        Inputs:
            kwargs:
                glyphRegistrationName: the name of the glyph, used to adjust the glyph outputs
        '''
        # parse input
        glyphRegistrationName = kwargs.get("glyphRegistrationName", "Glyph1")
        
        # get active view and source
        renderView1 = GetActiveViewOrCreate('RenderView')
        renderView1.UseColorPaletteForBackground = 0
        renderView1.Background = [1.0, 1.0, 1.0]

        # set source and field
        _source = "Transform1"
        _source_v = "Glyph1"
        field1 = "T"
        field2 = "viscosity"
        layout_resolution = (1350, 704)
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
       
        # find source
        source1 = FindSource(_source)
        sourceV = FindSource(_source_v)
    
        # show source1
        source1Display = Show(source1, renderView1, 'GeometryRepresentation')
        source1Display.SetScalarBarVisibility(renderView1, True)
        if PLOT_AXIS:
            source1Display.DataAxesGrid.GridAxesVisibility = 1
            source1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
        # get color transfer function/color map for 'field'
        field1LUT = GetColorTransferFunction(field1)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field1, 'Magnitude'))
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
        source1Display.SetScalarBarVisibility(renderView1, True)
        # Rescale transfer function, 2d transfer function
        field1LUT.RescaleTransferFunction(273.0, 2273.0)
        # colorbar position
        field1LUTColorBar = GetScalarBar(field1LUT, renderView1)
        field1LUTColorBar.Orientation = 'Horizontal'
        field1LUTColorBar.WindowLocation = 'Any Location'
        field1LUTColorBar.Position = [0.041, 0.908]
        field1LUTColorBar.ScalarBarLength = 0.33
        field1LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.TitleFontFamily = 'Times'
        field1LUTColorBar.LabelFontFamily = 'Times'
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        
        # show sourceV (vector field)
        sourceVDisplay = Show(sourceV, renderView1, 'GeometryRepresentation')
        sourceVDisplay.SetScalarBarVisibility(renderView1, True)
        # get color transfer function/color map for 'field'
        fieldVLUT = GetColorTransferFunction('velocity')
        if MAX_VELOCITY > 0.0:
            fieldVLUT.RescaleTransferFunction(0.0, MAX_VELOCITY)
        # colorbar position
        fieldVLUTColorBar = GetScalarBar(fieldVLUT, renderView1)
        fieldVLUTColorBar.Orientation = 'Horizontal'
        fieldVLUTColorBar.WindowLocation = 'Any Location'
        fieldVLUTColorBar.Position = [0.630, 0.908]
        fieldVLUTColorBar.ScalarBarLength = 0.33
        fieldVLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        fieldVLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        fieldVLUTColorBar.TitleFontFamily = 'Times'
        fieldVLUTColorBar.LabelFontFamily = 'Times'
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0


        # change point position
        pointName = "PointSource_" + glyphRegistrationName
        pointSource1 = FindSource(pointName)
        if "chunk" == "chunk":
            pointSource1.Center = [0, 6.7e6, 0]
        # show the representative point                                                                                                    
        _source_v_re = _source_v + "_representative"                                                                                       
        sourceVRE = FindSource(_source_v_re)                                                                                               
        sourceVREDisplay = Show(sourceVRE, renderView1, 'GeometryRepresentation') 
        # show the annotation
        _source_v_txt = _source_v + "_text" 
        sourceVTXT = FindSource(_source_v_txt)                                                                                               
        sourceVTXTDisplay = Show(sourceVTXT, renderView1, 'GeometryRepresentation')
        sourceVTXTDisplay.Color = [0.0, 0.0, 0.0]
        
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(layout_resolution[0], layout_resolution[1])
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [0.0, 5.6e5, 2.5e7]
            renderView1.CameraFocalPoint = [0.0, 6e6, 0.0]
            renderView1.CameraParallelScale = 8e5
        elif "GEOMETRY" == "box":
            renderView1.CameraPosition = [4700895.868280185, 2538916.5897593317, 15340954.822755022]
            renderView1.CameraFocalPoint = [4700895.868280185, 2538916.5897593317, 0.0]
            renderView1.CameraParallelScale = 487763.78047352127
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "T_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "T_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # second plot
        field2 = "viscosity"
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
        field2PWF = GetOpacityTransferFunction('viscosity')
        field2LUT.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        field2PWF.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field2, 'Magnitude'))
        source1Display.SetScalarBarVisibility(renderView1, True)
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        # Hide the scalar bar for the first field color map
        HideScalarBarIfNotNeeded(field1LUT, renderView1)
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(1350, 704)
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [0.0, 5.6e5, 2.5e7]
            renderView1.CameraFocalPoint = [0.0, 6e6, 0.0]
            renderView1.CameraParallelScale = 8e5
        elif "GEOMETRY" == "box":
            renderView1.CameraPosition = [4700895.868280185, 2538916.5897593317, 15340954.822755022]
            renderView1.CameraFocalPoint = [4700895.868280185, 2538916.5897593317, 0.0]
            renderView1.CameraParallelScale = 487763.78047352127
        # colorbar position
        field2LUTColorBar = GetScalarBar(field2LUT, renderView1)
        field2LUTColorBar.Orientation = 'Horizontal'
        field2LUTColorBar.WindowLocation = 'Any Location'
        field2LUTColorBar.Position = [0.041, 0.908]
        field2LUTColorBar.ScalarBarLength = 0.33
        field2LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.TitleFontFamily = 'Times'
        field2LUTColorBar.LabelFontFamily = 'Times'
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "viscosity_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "viscosity_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # hide plots
        Hide(source1, renderView1)
        Hide(sourceV, renderView1)
        Hide(sourceVRE, renderView1)
        Hide(sourceVTXT, renderView1)
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
        HideScalarBarIfNotNeeded(fieldVLUT, renderView1)
    
    def plot_step_whole(self, **kwargs): 
        '''
        plot a step
        Inputs:
            kwargs:
                glyphRegistrationName: the name of the glyph, used to adjust the glyph outputs
        '''
        # parse input
        glyphRegistrationName = kwargs.get("glyphRegistrationName", "Glyph1")
        
        # get active view and source
        renderView1 = GetActiveViewOrCreate('RenderView')
        renderView1.UseColorPaletteForBackground = 0
        renderView1.Background = [1.0, 1.0, 1.0]

        # set source and field
        _source = "Transform1"
        _source_v = "Glyph1"
        field1 = "T"
        field2 = "viscosity"
        layout_resolution = (1350, 704)
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
       
        # find source
        source1 = FindSource(_source)
        sourceV = FindSource(_source_v)
    
        # show source1
        source1Display = Show(source1, renderView1, 'GeometryRepresentation')
        source1Display.SetScalarBarVisibility(renderView1, True)
        if PLOT_AXIS:
            source1Display.DataAxesGrid.GridAxesVisibility = 1
            source1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
        # get color transfer function/color map for 'field'
        field1LUT = GetColorTransferFunction(field1)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field1, 'Magnitude'))
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
        source1Display.SetScalarBarVisibility(renderView1, True)
        # Rescale transfer function, 2d transfer function
        field1LUT.RescaleTransferFunction(273.0, 2273.0)
        # colorbar position
        field1LUTColorBar = GetScalarBar(field1LUT, renderView1)
        field1LUTColorBar.Orientation = 'Horizontal'
        field1LUTColorBar.WindowLocation = 'Any Location'
        field1LUTColorBar.Position = [0.041, 0.908]
        field1LUTColorBar.ScalarBarLength = 0.33
        field1LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.TitleFontFamily = 'Times'
        field1LUTColorBar.LabelFontFamily = 'Times'
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        
        # show sourceV (vector field)
        sourceVDisplay = Show(sourceV, renderView1, 'GeometryRepresentation')
        sourceVDisplay.SetScalarBarVisibility(renderView1, True)
        # get color transfer function/color map for 'field'
        fieldVLUT = GetColorTransferFunction('velocity')
        if MAX_VELOCITY > 0.0:
            fieldVLUT.RescaleTransferFunction(0.0, MAX_VELOCITY)
        # colorbar position
        fieldVLUTColorBar = GetScalarBar(fieldVLUT, renderView1)
        fieldVLUTColorBar.Orientation = 'Horizontal'
        fieldVLUTColorBar.WindowLocation = 'Any Location'
        fieldVLUTColorBar.Position = [0.630, 0.908]
        fieldVLUTColorBar.ScalarBarLength = 0.33
        fieldVLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        fieldVLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        fieldVLUTColorBar.TitleFontFamily = 'Times'
        fieldVLUTColorBar.LabelFontFamily = 'Times'
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0


        # change point position
        pointName = "PointSource_" + glyphRegistrationName
        pointSource1 = FindSource(pointName)
        if "chunk" == "chunk":
            pointSource1.Center = [0, 6.7e6, 0]
        # show the representative point                                                                                                    
        _source_v_re = _source_v + "_representative"                                                                                       
        sourceVRE = FindSource(_source_v_re)                                                                                               
        sourceVREDisplay = Show(sourceVRE, renderView1, 'GeometryRepresentation') 
        # show the annotation
        _source_v_txt = _source_v + "_text" 
        sourceVTXT = FindSource(_source_v_txt)                                                                                               
        sourceVTXTDisplay = Show(sourceVTXT, renderView1, 'GeometryRepresentation')
        sourceVTXTDisplay.Color = [0.0, 0.0, 0.0]
        
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(layout_resolution[0], layout_resolution[1])
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [1e5, 4.8e6, 2.5e7]
            renderView1.CameraFocalPoint = [1e5, 4.8e6, 0.0]
            renderView1.CameraParallelScale = 2.1e6
        elif "GEOMETRY" == "box":
            raise NotImplementError()
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "T_whole_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "T_whole_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # second plot
        field2 = "viscosity"
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
        field2PWF = GetOpacityTransferFunction('viscosity')
        field2LUT.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        field2PWF.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field2, 'Magnitude'))
        source1Display.SetScalarBarVisibility(renderView1, True)
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        # Hide the scalar bar for the first field color map
        HideScalarBarIfNotNeeded(field1LUT, renderView1)
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(1350, 704)
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [1e5, 4.8e6, 2.5e7]
            renderView1.CameraFocalPoint = [1e5, 4.8e6, 0.0]
            renderView1.CameraParallelScale = 2.1e6
        elif "GEOMETRY" == "box":
            raise NotImplementError()
        # colorbar position
        field2LUTColorBar = GetScalarBar(field2LUT, renderView1)
        field2LUTColorBar.Orientation = 'Horizontal'
        field2LUTColorBar.WindowLocation = 'Any Location'
        field2LUTColorBar.Position = [0.041, 0.908]
        field2LUTColorBar.ScalarBarLength = 0.33
        field2LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.TitleFontFamily = 'Times'
        field2LUTColorBar.LabelFontFamily = 'Times'
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "viscosity_whole_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "viscosity_whole_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # hide plots
        Hide(source1, renderView1)
        Hide(sourceV, renderView1)
        Hide(sourceVRE, renderView1)
        Hide(sourceVTXT, renderView1)
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
        HideScalarBarIfNotNeeded(fieldVLUT, renderView1)
    
    def plot_step_whole_whole(self, **kwargs): 
        '''
        plot a step
        Inputs:
            kwargs:
                glyphRegistrationName: the name of the glyph, used to adjust the glyph outputs
        '''
        # parse input
        glyphRegistrationName = kwargs.get("glyphRegistrationName", "Glyph1")
        
        # get active view and source
        renderView1 = GetActiveViewOrCreate('RenderView')
        renderView1.UseColorPaletteForBackground = 0
        renderView1.Background = [1.0, 1.0, 1.0]

        # set source and field
        _source = "Transform1"
        _source_v = "Glyph1"
        field1 = "T"
        field2 = "viscosity"
        layout_resolution = (1350, 704)
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
       
        # find source
        source1 = FindSource(_source)
    
        # show source1
        source1Display = Show(source1, renderView1, 'GeometryRepresentation')
        source1Display.SetScalarBarVisibility(renderView1, True)
        if PLOT_AXIS:
            source1Display.DataAxesGrid.GridAxesVisibility = 1
            source1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
            source1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
        # get color transfer function/color map for 'field'
        field1LUT = GetColorTransferFunction(field1)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field1, 'Magnitude'))
        HideScalarBarIfNotNeeded(field2LUT, renderView1)
        source1Display.SetScalarBarVisibility(renderView1, True)
        # Rescale transfer function, 2d transfer function
        field1LUT.RescaleTransferFunction(273.0, 2273.0)
        # colorbar position
        field1LUTColorBar = GetScalarBar(field1LUT, renderView1)
        field1LUTColorBar.Orientation = 'Horizontal'
        field1LUTColorBar.WindowLocation = 'Any Location'
        field1LUTColorBar.Position = [0.041, 0.908]
        field1LUTColorBar.ScalarBarLength = 0.33
        field1LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field1LUTColorBar.TitleFontFamily = 'Times'
        field1LUTColorBar.LabelFontFamily = 'Times'
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        
        # show sourceV (vector field), not included
        
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(layout_resolution[0], layout_resolution[1])
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [5e4, 2.8e6, 2.5e7]
            renderView1.CameraFocalPoint = [5e4, 2.8e6, 0.0]
            renderView1.CameraParallelScale = 5e6
        elif "GEOMETRY" == "box":
            raise NotImplementError()
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "T_whole_whole_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "T_whole_whole_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # second plot
        field2 = "viscosity"
        # get color transfer function/color map for 'field'
        field2LUT = GetColorTransferFunction(field2)
        field2PWF = GetOpacityTransferFunction('viscosity')
        field2LUT.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        field2PWF.RescaleTransferFunction(ETA_MIN, ETA_MAX)
        # set scalar coloring
        ColorBy(source1Display, ('POINTS', field2, 'Magnitude'))
        source1Display.SetScalarBarVisibility(renderView1, True)
        # hide the grid axis
        renderView1.OrientationAxesVisibility = 0
        # Hide the scalar bar for the first field color map
        HideScalarBarIfNotNeeded(field1LUT, renderView1)
        # adjust layout and camera & get layout & set layout/tab size in pixels
        layout1 = GetLayout()
        layout1.SetSize(1350, 704)
        renderView1.InteractionMode = '2D'
        if "GEOMETRY" == "chunk":
            renderView1.CameraPosition = [5e4, 2.8e6, 2.5e7]
            renderView1.CameraFocalPoint = [5e4, 2.8e6, 0.0]
            renderView1.CameraParallelScale = 5e6
        elif "GEOMETRY" == "box":
            raise NotImplementError()
        # colorbar position
        field2LUTColorBar = GetScalarBar(field2LUT, renderView1)
        field2LUTColorBar.Orientation = 'Horizontal'
        field2LUTColorBar.WindowLocation = 'Any Location'
        field2LUTColorBar.Position = [0.041, 0.908]
        field2LUTColorBar.ScalarBarLength = 0.33
        field2LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        field2LUTColorBar.TitleFontFamily = 'Times'
        field2LUTColorBar.LabelFontFamily = 'Times'
        # save figure
        fig_path = os.path.join(self.pv_output_dir, "viscosity_whole_whole_t%.4e.pdf" % self.time)
        fig_png_path = os.path.join(self.pv_output_dir, "viscosity_whole_whole_t%.4e.png" % self.time)
        SaveScreenshot(fig_png_path, renderView1, ImageResolution=layout_resolution)
        ExportView(fig_path, view=renderView1)

        # hide plots
        Hide(source1, renderView1)
        HideScalarBarIfNotNeeded(field2LUT, renderView1)


def add_glyph1(_source, field, scale_factor, **kwargs):
    '''
    add glyph in plots
    Inputs:
        scale_factor: scale of arrows
        ghost_field: the colorbar of a previous "ghost field" needs to be hide again to
            prevent it from being shown.
        kwargs:
            registrationName : the name of registration
            representative_value: a value to represent by the constant vector
    '''
    registrationName = kwargs.get("registrationName", 'Glyph1')
    representative_value = kwargs.get("representative_value", 0.05)
    
    # get active source and renderview
    pvd = FindSource(_source)
    renderView1 = GetActiveViewOrCreate('RenderView')

    # add glyph
    glyph1 = Glyph(registrationName=registrationName, Input=pvd, GlyphType='2D Glyph')
    # adjust orientation and scale
    glyph1.OrientationArray = ['POINTS', 'velocity']
    glyph1.ScaleArray = ['POINTS', 'velocity']
    glyph1.ScaleFactor = scale_factor
    glyph1.MaximumNumberOfSamplePoints = 20000

    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
    field0 = glyph1Display.ColorArrayName[1]
    field0LUT = GetColorTransferFunction(field0)
    # set the vector line width
    glyph1Display.LineWidth = 2.0
    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)
    # get color transfer function/color map for 'field'
    fieldLUT = GetColorTransferFunction(field)
    # get opacity transfer function/opacity map for 'field'
    fieldPWF = GetOpacityTransferFunction(field)

    # set scalar coloring
    ColorBy(glyph1Display, None)
    glyph1Display.AmbientColor = [1.0, 1.0, 1.0]
    glyph1Display.DiffuseColor = [1.0, 1.0, 1.0]

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(fieldLUT, renderView1)
    fieldLUT1 = GetColorTransferFunction(field0)
    HideScalarBarIfNotNeeded(fieldLUT1, renderView1)
    
    # add a representative vector
    pointName = "PointSource_" + registrationName
    pointSource1 = PointSource(registrationName=pointName)
    if "GEOMETRY" == "chunk":
        pointSource1.Center = [0, 6.4e6, 0]
    elif "GEOMETRY" == "box":
        pointSource1.Center = [4.65e6, 2.95e6, 0]
        
    # pointSource1.Center 
    pointSource1Display = Show(pointSource1, renderView1, 'GeometryRepresentation')
    print(dir(pointSource1))
    print(pointSource1.Center)

    calculatorName="Calculator_" + registrationName
    calculator1 = Calculator(registrationName=calculatorName, Input=pointSource1)
    calculator1.ResultArrayName = 'constant_velocity'
    calculator1.Function = '%.4e*iHat' % (scale_factor*representative_value)
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # add glyph
    glyph2Name = registrationName+"_representative"
    glyph2 = Glyph(registrationName=glyph2Name, Input=calculator1, GlyphType='2D Glyph')
    # adjust orientation and scale
    glyph2.OrientationArray = ['POINTS', 'constant_velocity']
    # glyph2.ScaleArray = ['POINTS', 'No scale array']
    glyph2.ScaleArray = ['POINTS', 'constant_velocity']
    # glyph2.ScaleFactor = 4e4
    glyph2.ScaleFactor = 1.0
    glyph2Display = Show(glyph2, renderView1, 'GeometryRepresentation')
    glyph2Display.AmbientColor = [0.0, 0.0, 0.0]
    glyph2Display.DiffuseColor = [0.0, 0.0, 0.0]
    # set the vector line width
    glyph2Display.LineWidth = 2.0
    # show color bar/color legend
    glyph2Display.SetScalarBarVisibility(renderView1, False)

    # add text
    # create a new 'Text'
    textName = registrationName + "_text"
    text1 = Text(registrationName=textName)
    # Properties modified on text1
    text1.Text = '5cm / yr'
    # show data in view
    text1Display = Show(text1, renderView1, 'TextSourceRepresentation')
    # Properties modified on text1Display
    text1Display.WindowLocation = 'Upper Center'
    text1Display.Color = [0.0, 0.0, 0.0]
    
    # hide data in view
    # Hide(pvd, renderView1)
    # hide glaph in view
    Hide(glyph1, renderView1)
    Hide(pointSource1, renderView1)
    Hide(calculator1)
    Hide(glyph2, renderView1)
    Hide(text1, renderView1)
    
    # update the view to ensure updated data information
    renderView1.Update()


def add_deformation_mechanism(_source, **kwargs):
    '''
    add programable filter to deal with the deformation mechanism
    Inputs:
        _source to extract the data
    kwargs:
        registrationName : the name of registration
    '''
    registrationName = kwargs.get("registrationName", 'Glyph1')
    
    # get active source and renderview
    pvd = FindSource(_source)
    renderView1 = GetActiveViewOrCreate('RenderView')

    # create a new 'Programmable Filter'
    programmableFilter1 = ProgrammableFilter(registrationName=registrationName, Input=pvd)
    programmableFilter1.Script = \
"""
import numpy as np
pressure = inputs[0].PointData["p"].GetArrays()[0]
strain_rate = inputs[0].PointData["strain_rate"].GetArrays()[0]
diff_visc = inputs[0].PointData["diffusion_viscosity"].GetArrays()[0]
disl_visc= inputs[0].PointData["dislocation_viscosity"].GetArrays()[0]
p_visc= inputs[0].PointData["peierls_viscosity"].GetArrays()[0]
y_visc1 = (pressure * np.sin(25.0/180.0 * np.pi) + 50e6) / 2.0 / strain_rate
y_visc2 = 500e6 / 2.0 / strain_rate
y_visc = np.minimum(y_visc1, y_visc2)
d_mech = np.argmin(np.array([diff_visc, disl_visc, y_visc, p_visc]), axis=0)
# output.PointData.append(d_mech, "deformation_mechanism")
"""
    programmableFilter1.RequestInformationScript = ''
    programmableFilter1.RequestUpdateExtentScript = ''
    programmableFilter1.PythonPath = ''



def main():
    all_available_graphical_snapshots = ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS
    all_available_graphical_times = ALL_AVAILABLE_GRAPHICAL_TIMES
    # types of plots included
    plot_types = PLOT_TYPES
    assert(len(all_available_graphical_snapshots) == len(all_available_graphical_times))
    # First, make directory for images if it's not there
    if not os.path.isdir("IMG_OUTPUT_DIR"):
        os.mkdir("IMG_OUTPUT_DIR")
    # Process and generate plots
    Slab = SLAB("PARAVIEW_FILE", output_dir="IMG_OUTPUT_DIR")
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
                if "upper_mantle" in plot_types:
                    Slab.plot_step()
                if "whole" in plot_types:
                    Slab.plot_step_whole()
                if "whole_whole" in plot_types:
                    Slab.plot_step_whole_whole()
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
