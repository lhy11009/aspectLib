# variable to substitute
# vts file:
#   PARAVIEW_FILE
# directory for images:
#   IMG_OUTPUT_DIR

class SLAB(PARAVIEW_PLOT):
    '''
    Inherit frome PARAVIEW_PLOT
    Usage is for plotting the shape of the 3D slabs
    '''
    def plot_surface_slice(self):
        '''
        Generate a visualization of a slice perpendicular to z direction, and locates
        at 10 km in depth.
        '''
        slice1, slice1Display, _ = add_slice(self.solutionpvd, "sp_upper", [2000000.0, 2000000.0, 990000.0],\
        [0.0, 0.0, 1.0], renderView=self.renderView1)
        adjust_camera(self.renderView1, [2000000.0, 2000000.0, 14540803.753676033],\
        [2000000.0, 2000000.0, 990000.0], 4019667.7972010756, [0.0, 1.0, 0.0])
        file_out = os.path.join(self.output_dir, "surface_slice.png")
        SaveScreenshot(file_out, self.renderView1, ImageResolution=[1148, 792])


def main():
    # First, make directory for images if it's not there
    if not os.path.isdir("IMG_OUTPUT_DIR"):
        os.mkdir("IMG_OUTPUT_DIR")
    # Process and generate plots
    Slab = SLAB("PARAVIEW_FILE", output_dir="IMG_OUTPUT_DIR")
    Slab.plot_surface_slice()


main()
