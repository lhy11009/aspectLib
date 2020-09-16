# variable to substitute
# VISIT_PARTICLE_FILE:
#   vts file
# DATA_OUTPUT_DIR
#   directory to output data
OpenDatabase("VISIT_PARTICLE_FILE", 0)
AddPlot("Mesh", "mesh", 1, 1)
AddPlot("Molecule", "id", 1, 1)
DrawPlots()
ExportDBAtts = ExportDBAttributes()
ExportDBAtts.allTimes = 1
ExportDBAtts.dirname = "DATA_OUTPUT_DIR"
ExportDBAtts.filename = "visit_particles"
ExportDBAtts.timeStateFormat = "_%04d"
ExportDBAtts.db_type = "XYZ"
ExportDBAtts.db_type_fullname = "XYZ_1.0"
ExportDBAtts.variables = ("id")
ExportDBAtts.writeUsingGroups = 0
ExportDBAtts.groupSize = 48
ExportDatabase(ExportDBAtts)
