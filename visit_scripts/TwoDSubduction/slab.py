# variable to substitute
# vts file:
#   VISIT_FILE
# initial adaptive refinement level:
#   INITIAL_ADAPTIVE_REFINEMENT
# directory for images:
#   IMG_OUTPUT_DIR
# directory for output particle information
#   PARTICLE_OUTPUT_DIR
# all available snapshots for graphical output
#   ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS
# all available snapshots for particle output
#   ALL_AVAILABLE_PARTICLE_SNAPSHOTS


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
        self.filein = filein
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
    
    def set_view_attrs(self, window_coords):
        '''
        set the attributes of camera
        '''
        # assertion
        assert len(window_coords) == 4
        # fix camera
        # Begin spontaneous state
        View2DAtts = View2DAttributes()
        View2DAtts.windowCoords = window_coords
        View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
        View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
        View2DAtts.fullFrameAutoThreshold = 100
        View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.windowValid = 1
        SetView2D(View2DAtts)
        # End spontaneous state
    
    def set_rotation_transform(self, rot_deg):
        '''
        set the rotation
        '''
        AddOperator("Transform", 1)
        TransformAtts = TransformAttributes()
        TransformAtts.doRotate = 1
        TransformAtts.rotateOrigin = (0, 0, 0)
        TransformAtts.rotateAxis = (0, 0, 1)
        TransformAtts.rotateAmount = rot_deg
        TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
        TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
        TransformAtts.outputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
        SetOperatorOptions(TransformAtts, 0, 1)
        pass

    def set_threshold(self, name_, values):
        '''
        set values of threshhold
        Inputs:
            name_(str): name of the field
            values(list of 2): lower and upper limit
        '''
        SetActivePlots((self.idxs[name_]))
        ThresholdAtts = ThresholdAttributes()
        ThresholdAtts.outputMeshType = 0
        ThresholdAtts.boundsInputType = 0
        ThresholdAtts.listedVarNames = ("default")
        ThresholdAtts.zonePortions = (1)
        ThresholdAtts.lowerBounds = (values[0])
        ThresholdAtts.upperBounds = (values[1])
        ThresholdAtts.defaultVarName = name_
        ThresholdAtts.defaultVarIsScalar = 1
        ThresholdAtts.boundsRange = ("%.6e:%.6e" % (values[0], values[1]))
        SetOperatorOptions(ThresholdAtts, 1, 0)
        DrawPlots()
    
    def set_pseudo_color(self, name_, **kwargs):
        '''
        set pseudocolor plot
        '''
        SetActivePlots(self.idxs[name_])
        PseudocolorAtts = PseudocolorAttributes()
        PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot

        # color
        color_table = kwargs.get('color_table', "viridis")
        PseudocolorAtts.colorTableName = color_table

        # min and max values to use
        try:
            limits = kwargs['limits']
        except KeyError:
            pass
        else:
            assert len(limits) == 2
            PseudocolorAtts.minFlag = 1
            PseudocolorAtts.min = limits[0]
            PseudocolorAtts.useBelowMinColor = 0
            PseudocolorAtts.belowMinColor = (0, 0, 0, 255)
            PseudocolorAtts.maxFlag = 1
            PseudocolorAtts.max = limits[1]
            PseudocolorAtts.useAboveMaxColor = 0
            PseudocolorAtts.aboveMaxColor = (0, 0, 0, 255)
        
        # invert color table
        invert_color = kwargs.get('invert_color', False)
        if invert_color: 
            PseudocolorAtts.invertColorTable = 1

        # use log value
        log = kwargs.get('log', False)
        if log:
            PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
        
        # commit settings 
        SetPlotOptions(PseudocolorAtts)
        
    def save_window(self, filename, **kwargs):
        '''
        save window
        Inputs:)
            filename(str): base name for file
            kwargs(dict):
                size(list of 2): size of plot in (width, height)
        '''
        SaveWindowAtts = SaveWindowAttributes()
        SaveWindowAtts.outputDirectory = self.output_dir
        SaveWindowAtts.family = 0  # we don't use family name here
        SaveWindowAtts.fileName = "%s/%s%06d" % (self.output_dir, filename, self.time_snap)
        SaveWindowAtts.format = SaveWindowAtts.PNG

        # size
        try:
            size = kwargs['size']
        except KeyError:
            pass
        else:
            assert len(size) == 2
            SaveWindowAtts.width = size[0]
            SaveWindowAtts.height = size[1]
        
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()
    
    def abort(self):
        '''
        abort mission and leave
        '''
        # delete active plots
        all_idxs_tuple = tuple([i for i in range(self.all_idxs)])
        SetActivePlots(all_idxs_tuple)
        DeleteActivePlots()

        # close file 
        OpenDatabase(self.filein)
    

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
        plot_types, vars_ = self.get_plots()

        # call __init__ function of parent
        VISIT_PLOT.__init__(self, filein, plot_types=plot_types, vars=vars_, output_dir=output_dir)

        # define the field of the slab
        DefineScalarExpression("slab", "spcrust+spharz")
        
        # define a new field if needed
        if IF_DEFORM_MECHANISM:
            DefineScalarExpression("deform_mechanism", "if(lt(viscosity, 0.99e24),if(lt(dislocation_viscosity,diffusion_viscosity), 1.0, 0.0), 2.0)")
        self.add_plot("Pseudocolor", "deform_mechanism")
        
        # set transformation
        self.set_rotation_transform(52.0)

        # set thresholds
        threshold_idxs_tuple=tuple([self.idxs['spcrust'], self.idxs['spharz']])
        SetActivePlots(threshold_idxs_tuple)
        AddOperator("Threshold", 0)
    
        # draw plots
        all_idxs_tuple = tuple([i for i in range(self.all_idxs)])
        SetActivePlots(all_idxs_tuple)
        DrawPlots()
        HideActivePlots()
    
    def get_plots(self):
        """
        get types of plot
        to be reloaded in children
        Returns:
            plot_types(list)
            vars_(list): variables to plot
        """
        plot_types = ["Mesh", "Pseudocolor", "Pseudocolor", "Pseudocolor", "Pseudocolor", "Vector"]
        vars_ = ["mesh", "spcrust", "spharz", "T", "viscosity", "velocity"]
        return plot_types, vars_

    def plot_time_snap(self):
        """
        plot a single time step
        """
        # spcrust
        self.plot_crust()
        
        # plot slab viscosity
        self.plot_viscosity_slab()

        # plot deform_mechanism 
        if IF_DEFORM_MECHANISM:
            self.plot_deform_mechanism_slab()

        # plot upper mantle
        # viscosity
        self.plot_viscosity_upper_mantle()

        # deform mechanism 
        if IF_DEFORM_MECHANISM:
            self.plot_deform_mechanism_upper_mantle()

    
    def plot_crust(self):
        '''
        plot crust properties
        '''
        # set camera
        self.set_view_attrs((-200000, 200000, 6.1e+06, 6.372e+06))

        # threshold for spcrust
        self.set_threshold('spcrust', [0.8, 1e+37])

        # threshold for spharz
        self.set_threshold('spharz', [0.8, 1e+37])
        
        # set min and max value
        self.set_pseudo_color('spcrust', color_table='Reds', limits=[0.0, 1.0])
        self.set_pseudo_color('spharz', color_table='Blues', limits=[0.0, 1.0])

        # show active plots
        SetActivePlots((self.idxs['mesh'], self.idxs['spcrust'], self.idxs['spharz']))
        HideActivePlots()
        
        # Save Plot
        self.save_window('slab_crust_snap', size=[2048, 2048])
        HideActivePlots()

    def plot_viscosity_slab(self):
        '''
        plot viscosity
        '''
        # set camera
        self.set_view_attrs((-200000, 200000, 6.1e+06, 6.372e+06))
        
        # set up viscosity
        # change to log scale and invert the color table
        self.set_pseudo_color('viscosity', color_table="SCM_romaO", invert_color=True, log=True)
       
        # set up velocity
        SetActivePlots(self.idxs['velocity'])
        VectorAtts = VectorAttributes()
        VectorAtts.glyphLocation = VectorAtts.AdaptsToMeshResolution  # AdaptsToMeshResolution, UniformInSpace
        VectorAtts.colorTableName = "BrBG"
        VectorAtts.nVectors = 4000
        VectorAtts.lineWidth = 0
        VectorAtts.scale = 0.025
        VectorAtts.scaleByMagnitude = 1
        VectorAtts.autoScale = 1
        SetPlotOptions(VectorAtts)
        
        SetActivePlots((self.idxs['viscosity'], self.idxs['velocity']))
        HideActivePlots()

        # save plot 
        self.save_window('slab_viscosity_snap')
        HideActivePlots()
    
    def plot_deform_mechanism_slab(self):
        '''
        plot deform mechanism
        '''
        # set camera
        self.set_view_attrs((-200000, 200000, 6.1e+06, 6.372e+06))
        
        # set up deform_mechanism
        self.set_pseudo_color('deform_mechanism', color_table='viridis')
       
        SetActivePlots((self.idxs['deform_mechanism']))
        HideActivePlots()

        # save plot 
        self.save_window('slab_deform_mechanism_snap')
        HideActivePlots()

    def plot_viscosity_upper_mantle(self):
        '''
        plot viscosity
        '''

        # set camera
        self.set_view_attrs((-1.0e+06, 1.0e+06, 5.4e+06, 6.4e+06))
        
        # set up viscosity
        # change to log scale and invert the color table
        self.set_pseudo_color('viscosity', color_table="SCM_romaO", invert_color=True, log=True)
       
        # set up velocity
        SetActivePlots(self.idxs['velocity'])
        VectorAtts = VectorAttributes()
        VectorAtts = VectorAttributes( )
        VectorAtts.glyphLocation = VectorAtts.UniformInSpace # AdaptsToMeshResolution, UniformInSpace
        VectorAtts.colorTableName = "BrBG"
        VectorAtts.nVectors = 20000
        VectorAtts.lineWidth = 0
        VectorAtts.scale = 0.05
        VectorAtts.scaleByMagnitude = 1
        VectorAtts.autoScale = 1
        SetPlotOptions(VectorAtts)

        SetActivePlots((self.idxs['viscosity'], self.idxs['velocity']))
        HideActivePlots()
        # save plot 
        self.save_window('um_viscosity_snap')
        HideActivePlots()
    
    def plot_deform_mechanism_upper_mantle(self):
        '''
        plot deform mechanism in upper mantle
        '''
        # set camera
        self.set_view_attrs((-1.0e+06, 1.0e+06, 5.4e+06, 6.4e+06))
        
        # set up deform_mechanism
        self.set_pseudo_color('deform_mechanism', color_table='viridis')
       
        SetActivePlots((self.idxs['deform_mechanism']))
        HideActivePlots()

        # save plot 
        self.save_window('um_deform_mechanism_snap')
        HideActivePlots()
        

class EXPORT_PARTICLE(VISIT_PLOT):
    '''
    export particles from visit file
    generate a file with positions of particles for each step
    '''
    
    def __init__(self, filein, **kwargs):
        """
        initiate
        Args:
            filein(str): path of vtu file
        """
        # get options
        output_dir = kwargs.get('output_dir', '.')

        # types of plot
        plot_types, vars_ = self.get_plots()

        # call __init__ function of parent
        VISIT_PLOT.__init__(self, filein, plot_types=plot_types, vars=vars_, output_dir=output_dir)

        # draw plots
        all_idxs_tuple = tuple([i for i in range(self.all_idxs)])
        SetActivePlots(all_idxs_tuple)
        DrawPlots()
        HideActivePlots()
    
    def get_plots(self):
        """
        get types of plot
        to be reloaded in children
        Returns:
            plot_types(list)
            vars_(list): variables to plot
        """
        plot_types = ["Mesh", "Molecule"]
        vars_ = ["mesh", "id"]
        return plot_types, vars_
    
    def plot_time_snap(self):
        # set active
        SetActivePlots((self.idxs['id']))
        # SetActivePlots((self.idxs['mesh'], self.idxs['id']))
        HideActivePlots()

        # export 
        ExportDBAtts = ExportDBAttributes()
        ExportDBAtts.allTimes = 1
        ExportDBAtts.dirname = "PARTICLE_OUTPUT_DIR"
        ExportDBAtts.filename = "visit_particles_%06d" % self.time_snap
        ExportDBAtts.timeStateFormat = "_%04d"
        ExportDBAtts.db_type = "XYZ"
        ExportDBAtts.db_type_fullname = "XYZ_1.0"
        ExportDBAtts.variables = ("id")
        ExportDBAtts.writeUsingGroups = 0
        ExportDBAtts.groupSize = 48
        ExportDatabase(ExportDBAtts)
        HideActivePlots()



def main():
    if IF_PLOT_SLAB:
        Slab = SLAB("VISIT_FILE", output_dir="IMG_OUTPUT_DIR")
        # First number is the number of initial adaptive refinements
        # Second one is the snapshot to plot
        # here we prefer to use a series of snapshots.
        # If this doesn't work, we will use a single snapshot
        steps = PLOT_SLAB_STEPS
        if not steps == []:
            for step in steps:
                Slab(INITIAL_ADAPTIVE_REFINEMENT+step)
        else:
            Slab(INITIAL_ADAPTIVE_REFINEMENT+SINGLE_SNAPSHOT)
        Slab.abort()

    # Output particles for slab morphology
    if IF_EXPORT_SLAB_MORPH:
        Export_Particle = EXPORT_PARTICLE("VISIT_PARTICLE_FILE", output_dir="DATA_OUTPUT_DIR")
        for snapshot in ALL_AVAILABLE_PARTICLE_SNAPSHOTS:
            Export_Particle(snapshot)
        Export_Particle.abort()

main()