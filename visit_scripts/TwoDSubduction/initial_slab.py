# variable to substitute
# VISIT_FILE:
#   vts file
# INITIAL_ADAPTIVE_REFINEMENT:
#   initial adaptive refinement level
# IMG_OUTPUT_DIR
#   directory for images
OpenDatabase("localhost:VISIT_FILE", 0)
AddPlot("Mesh", "mesh", 1, 1)
AddPlot("Pseudocolor", "spcrust", 1, 1)
DrawPlots()
for _ in range(INITIAL_ADAPTIVE_REFINEMENT):
    TimeSliderNextState()
DrawPlots()
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
# Same Plot
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputDirectory = "IMG_OUTPUT_DIR"
SaveWindowAtts.fileName = "IMG_OUTPUT_DIR/visit_initial_slab"
SaveWindowAtts.format = SaveWindowAtts.PNG
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
