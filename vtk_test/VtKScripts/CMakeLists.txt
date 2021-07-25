
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(VtKScripts)

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
  QUIET
)

if (NOT VTK_FOUND)
  message("Skipping VtKScripts: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(VtKScripts MACOSX_BUNDLE VtKScripts.cxx )
  target_link_libraries(VtKScripts PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(VtKScripts MACOSX_BUNDLE VtKScripts.cxx )
  target_link_libraries(VtKScripts PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS VtKScripts
    MODULES ${VTK_LIBRARIES}
    )
endif ()
