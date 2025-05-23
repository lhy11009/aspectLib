
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(AspectVTK)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkIOXML
  vtkIOLegacy
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkCommonDataModel
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersModeling
  vtkFiltersGeometry
  QUIET
)

include_directories("${PROJECT_BINARY_DIR}")

if (NOT VTK_FOUND)
  message("Skipping AspectVTK: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(AspectVTK MACOSX_BUNDLE AspectVTK.cxx )
  target_link_libraries(AspectVTK PRIVATE ${VTK_LIBRARIES})
else ()
  # build SlabAnalysis
  add_executable(TwoDSubduction_SlabAnalysis MACOSX_BUNDLE TwoDSubduction_SlabAnalysis.cxx AspectVTK.cxx )
  target_link_libraries(TwoDSubduction_SlabAnalysis PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS TwoDSubduction_SlabAnalysis
    MODULES ${VTK_LIBRARIES}
    )
  # build MOW
  add_executable(TwoDSubduction_MOW MACOSX_BUNDLE TwoDSubduction_MOW.cxx AspectVTK.cxx )
  target_link_libraries(TwoDSubduction_MOW PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS TwoDSubduction_MOW
    MODULES ${VTK_LIBRARIES}
    )
endif ()
