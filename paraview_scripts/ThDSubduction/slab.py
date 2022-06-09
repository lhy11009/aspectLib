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
    def setup_surface_slice(self):
        '''
        Generate a visualization of a slice perpendicular to z direction, and locates
        at 10 km in depth.
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, 2000000.0, 990000.0],\
        [0.0, 0.0, 1.0], renderView=self.renderView1, name="surface_z")
        adjust_camera(self.renderView1, [2000000.0, 2000000.0, 14540803.753676033],\
        [2000000.0, 2000000.0, 990000.0], 4019667.7972010756, [0.0, 1.0, 0.0])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view
    
    def setup_trench_slice_center(self):
        '''
        Generate a visualization of a slice perpendicular to y direction, and locates at trench center (y=0)
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, 1.0, 500000.0],\
        [0.0, 1.0, 0.0], renderView=self.renderView1, name="trench_center_y")
        adjust_camera(self.renderView1, [2000000.0, 14540803.753676033, 500000.0],\
        [2000000.0, 1.0, 500000.0], 4019667.7972010756, [0.0, 0.0, 1.0])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view

    def setup_trench_slice_edge(self):
        '''
        Generate a visualization of a slice perpendicular to y direction, and locates at trench center (y=0)
        Here I first set the plots up. This is equivalent to add filters to a paraview window in the gui.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, TRENCH_EDGE_Y, 500000.0],\
        [0.0, 1.0, 0.0], renderView=self.renderView1, name="trench_edge_y")
        adjust_camera(self.renderView1, [2000000.0, 14540803.753676033, 500000.0],\
        [2000000.0, TRENCH_EDGE_Y, 500000.0], 4019667.7972010756, [0.0, 0.0, 1.0])
        Hide3DWidgets()  # this is the same thing as unchecking the "show plane"
        Hide(slice1, self.renderView1) # hide data in view
    
    def setup_slab_iso_volume_upper(self):
        '''
        Generate a visualization of the iso volume of the upper crust
        '''
        isoVolume1, isoVolume1Display, _ = add_isovolume(self.solutionpvd, "sp_upper", (0.8, 1.0), name="slab_upper")
        Hide(isoVolume1, self.renderView1)  # hide data in view

 
    def plot_slice(self, filename_base):
        '''
        Plot surface slice
        After the the plots are setup, this function is then called to execute the exportation of figures.
        '''
        field_name = 'sp_upper'
        slice1 = FindSource(field_name)
        SetActiveSource(slice1)
        renderView1 = GetActiveViewOrCreate('RenderView') 
        # Show(slice1, renderView1, 'GeometryRepresentation')
        # adjust colorbar and camera
        sp_upperLUT = GetColorTransferFunction(field_name)
        adjust_slice_colorbar_camera(self.renderView1, sp_upperLUT)
        # save figure
        file_out = os.path.join(self.output_dir, "%s_%.4e.png" % (filename_base, self.time))
        SaveScreenshot(file_out, self.renderView1, ImageResolution=[1148, 792])
        print("Figure saved: %s" % file_out)


def adjust_slice_colorbar_camera(renderView, colorLUT):
    '''
    adjust colorbar and camera for a slice
    Inputs:
        renderView: an instance of the rendered view
        colorLUT: an instance of the colorbar
    '''
    # adjust colorbar
    colorLUTColorBar = GetScalarBar(colorLUT, renderView)
    colorLUTColorBar.WindowLocation = 'Any Location'
    colorLUTColorBar.ScalarBarLength = 0.33000000000000007
    colorLUTColorBar.Orientation = 'Horizontal'
    colorLUTColorBar.Position = [0.3458232931726907, 0.2540226986128623]
    # adjust camera
    adjust_camera(renderView, [2000000.0, 14540803.753676033, 500000.0],\
    [2000000.0, 1.0, 500000.0], 1875204.6933857028, [0.0, 0.0, 1.0])


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
                Slab.plot_slice("trench_slice_edge")
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
