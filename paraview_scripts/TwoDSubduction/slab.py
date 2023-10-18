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
        add_glyph1("Transform1", "velocity", "T", 5e6, registrationName="Glyph1")

    def plot_step(self): 
        '''
        plot a step
        '''
        pass

def main():
    all_available_graphical_snapshots = ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS
    all_available_graphical_times = ALL_AVAILABLE_GRAPHICAL_TIMES
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
