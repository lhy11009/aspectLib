# variable to substitute
# vts file:
#   VISIT_FILE
# initial adaptive refinement level:
#   INITIAL_ADAPTIVE_REFINEMENT
# directory for images:
#   IMG_OUTPUT_DIR


class VISIT_PLOT():
    """
    BASE CLASS for visit plot()
    Attributes:
        idx(dict): a dictionary for index
        all_idxs(int): total indexes
        time_snap(int): current time snap number
        output_dir(str): path of output
    """

    def __init__(self, filein, **kwargs):
        """
        initiate
        Args:
            filein(str): path of vtu file
            plots(list): lists of plot
        """
        # open data base
        OpenDatabase("localhost:%s" % filein, 0)

        # assign initial values 
        self.idxs = {}
        self.all_idxs = 0
        self.time_snap = 0
        self.output_dir = kwargs.get('output_dir', '.')

        # add plots
        # self.set_plots() is to be reloaded in children
        plot_types, vars_ = self.get_plots()
        for i in range(len(plot_types)):
            plot_type = plot_types[i]
            var_ = vars_[i]
            self.add_plot(plot_type, var_)

    def add_plot(self, plot_type, var):
        """
        add a plot
        Args:
            plot_type(str): type of plot
            var(str): variable to plot
        """
        AddPlot(plot_type, var, 1, 1)
        self.idxs[var] = self.all_idxs
        self.all_idxs += 1

    def goto_time_snap(self, time_snap):
        """
        goto time_snap
        Args:
            time_snap(int): number of time snap
        """
        for _ in range(time_snap - self.time_snap):
            TimeSliderNextState()
        self.time_snap = time_snap

    def __call__(self, time_snap):
        """
        call
        Args:
            time_snap(int): time_snap to plot
        """
        # shift to specific time snap
        self.goto_time_snap(time_snap)
        self.plot_time_snap()

    def plot_time_snap(self):
        """
        plot a single time step
        to be reload in children
        """
        pass

    def get_plots(self):
        """
        get types of plot
        to be reloaded in children
        Returns:
            plot_types(list)
            vars_(list): variables to plot
        """
        plot_types = []
        vars_ = []
        return plot_types, vars_


class SLAB(VISIT_PLOT):

    def __init__(self, filein, **kwargs):
        """
        initiate
        Args:
            filein(str): path of vtu file
        """
        # get options
        output_dir = kwargs.get('output_dir', '.')

        # types of plot
        plot_types = ["Mesh", "Pseudocolor", "Pseudocolor", "Pseudocolor", "Pseudocolor"]
        vars_ = ["mesh", "spcrust", "spharz", "T", "viscosity"]

        # call __init__ function of parent
        VISIT_PLOT.__init__(self, filein, plot_types=plot_types, vars=vars_, output_dir=output_dir)
    
    def get_plots(self):
        """
        get types of plot
        to be reloaded in children
        Returns:
            plot_types(list)
            vars_(list): variables to plot
        """
        plot_types = ["Mesh", "Pseudocolor", "Pseudocolor", "Pseudocolor", "Pseudocolor"]
        vars_ = ["mesh", "spcrust", "spharz", "T", "viscosity"]
        return plot_types, vars_

    def plot_time_snap(self):
        """
        plot a single time step
        """
        # fix camera
        # Begin spontaneous state
        View2DAtts = View2DAttributes()
        View2DAtts.windowCoords = (4.75068e+06, 5.01824e+06, 3.76314e+06, 4.02705e+06)
        View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
        View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
        View2DAtts.fullFrameAutoThreshold = 100
        View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.windowValid = 1
        SetView2D(View2DAtts)
        # End spontaneous state
        
        # draw plots
        all_idxs_tuple = tuple([i for i in range(self.all_idxs)])
        SetActivePlots(all_idxs_tuple)
        DrawPlots()
        
        # spcrust
        HideActivePlots()
        SetActivePlots((self.idxs['mesh'], self.idxs['spcrust']))
        HideActivePlots()
        # set min and max value
        PseudocolorAtts = PseudocolorAttributes()
        PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
        PseudocolorAtts.minFlag = 1
        PseudocolorAtts.min = 0
        PseudocolorAtts.useBelowMinColor = 0
        PseudocolorAtts.belowMinColor = (0, 0, 0, 255)
        PseudocolorAtts.maxFlag = 1
        PseudocolorAtts.max = 1
        PseudocolorAtts.useAboveMaxColor = 0
        PseudocolorAtts.aboveMaxColor = (0, 0, 0, 255)
        SetPlotOptions(PseudocolorAtts)
        
        # Save Plot
        SaveWindowAtts = SaveWindowAttributes()
        SaveWindowAtts.outputDirectory = self.output_dir
        SaveWindowAtts.fileName = "%s/visit_initial_slab_crust" % self.output_dir
        SaveWindowAtts.format = SaveWindowAtts.PNG
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()
        
        # Plot viscosity
        HideActivePlots()
        SetActivePlots((self.idxs['viscosity']))
        HideActivePlots()
        
        # change to log scale
        PseudocolorAtts = PseudocolorAttributes()
        PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
        SetPlotOptions(PseudocolorAtts)
        SaveWindowAtts.outputDirectory = self.output_dir
        SaveWindowAtts.fileName = "%s/visit_initial_slab_viscosity" % self.output_dir
        SaveWindowAtts.format = SaveWindowAtts.PNG
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()


def main():
    Slab = SLAB("VISIT_FILE", output_dir="IMG_OUTPUT_DIR")
    Slab(INITIAL_ADAPTIVE_REFINEMENT)


main()