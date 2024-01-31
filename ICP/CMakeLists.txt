set(classes
  vtkPolyDataCorrespondenceFilter
  vtkICP
  vtkMultiMeshAlignment
  vtkImplicitPolyDataDistance2
)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/icp/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/icp/config.h" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

set(private_headers
  icp/icp.h
  icp/point2point.h
  icp/point2plane.h
  icp/plane2plane.h
  icp/kabsch.h
)

# If we add classes, we will get an install warning, the subdirectorie icp
# is not installed. We add the source instead as pure source files
# This is new in CMake 3.23
set(nowrap_classes
  icp/icp
  icp/point2plane
  icp/point2point
  icp/plane2plane
  icp/kabsch)

# Source files (instead of classes)
set(sources
  icp/icp.cxx
  icp/point2plane.cxx
  icp/point2point.cxx
  icp/plane2plane.cxx
  icp/kabsch.cxx
  multiMeshInternals.cxx)
  
if(DEFINED VTK_MODULE_ENABLE_VTK_ICP)
  vtk_module_add_module(VTK::ICP
    CLASSES           ${classes}
    TEMPLATE_CLASSES  ${template_classes}
    SOURCES           ${sources}
    #NOWRAP_CLASSES    ${nowrap_classes} # Gives a warning in CMake 3.23
    HEADERS           ${headers}
    PRIVATE_TEMPLATES ${private_templates}
    #PRIVATE_HEADERS   ${private_headers} # Gives a warning in CMake 3.23
  )
  
  vtk_module_link(VTK::ICP
    PRIVATE
    VTK::CommonDataModel)
  
endif()
  
