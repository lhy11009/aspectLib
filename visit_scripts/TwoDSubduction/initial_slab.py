# variable to substitute
# vts file:
#   VISIT_FILE
# initial adaptive refinement level:
#   INITIAL_ADAPTIVE_REFINEMENT
# directory for images:
#   IMG_OUTPUT_DIR
OpenDatabase("localhost:VISIT_FILE", 0)

# add plots
# index of each plot
idxs={}
AddPlot("Mesh", "mesh", 1, 1)
idxs['mesh'] = 0
AddPlot("Pseudocolor", "spcrust", 1, 1)
idxs['spcrust'] = 1
AddPlot("Pseudocolor", "viscosity", 1, 1)
idxs['viscosity'] = 2
all_idxs = tuple([ item[1] for item in idxs.items() ])

# shift to the initial timestep
for _ in range(INITIAL_ADAPTIVE_REFINEMENT):
    TimeSliderNextState()

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
SetActivePlots(all_idxs)
DrawPlots()

# spcrust
HideActivePlots()
SetActivePlots((idxs['mesh'], idxs['spcrust']))
HideActivePlots()

# Save Plot
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputDirectory = "IMG_OUTPUT_DIR"
SaveWindowAtts.fileName = "IMG_OUTPUT_DIR/visit_initial_slab_crust"
SaveWindowAtts.format = SaveWindowAtts.PNG
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

# Plot viscosity
HideActivePlots()
SetActivePlots((idxs['viscosity']))
HideActivePlots()

# change to log scale
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
SetPlotOptions(PseudocolorAtts)
SaveWindowAtts.outputDirectory = "IMG_OUTPUT_DIR"
SaveWindowAtts.fileName = "IMG_OUTPUT_DIR/visit_initial_slab_viscosity"
SaveWindowAtts.format = SaveWindowAtts.PNG
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
