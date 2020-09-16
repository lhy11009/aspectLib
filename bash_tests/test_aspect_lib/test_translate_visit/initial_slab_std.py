# variable to substitute
# /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12/output1/solution.visit:
#   vts file
# 6:
#   initial adaptive refinement level
# /home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12/img
#   directory for images
OpenDatabase("localhost:/home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12/output1/solution.visit", 0)
AddPlot("Mesh", "mesh", 1, 1)
AddPlot("Pseudocolor", "spcrust", 1, 1)
DrawPlots()
for _ in range(6):
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
SaveWindowAtts.outputDirectory = "/home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12/img"
SaveWindowAtts.fileName = "/home/lochy/ASPECT_PROJECT/TwoDSubduction/isosurf_global2/isosurfULV3.000e+01testS12/img/visit_initial_slab"
SaveWindowAtts.format = SaveWindowAtts.PNG
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()