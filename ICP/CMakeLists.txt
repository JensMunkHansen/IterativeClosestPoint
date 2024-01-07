set(classes
  vtkPolyDataCorrespondenceFilter
  vtkICP
  vtkImplicitPolyDataDistance2
)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/icp/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/icp/config.h" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

set(private_headers
  icp/icp.h
  icp/point2point.h
  icp/point2plane.h
  icp/plane2plane.h
)
set(nowrap_classes
  icp/icp
  icp/point2plane
  icp/point2point
  icp/plane2plane)

if(DEFINED VTK_MODULE_ENABLE_VTK_ICP)
  vtk_module_add_module(VTK::ICP
    CLASSES           ${classes}
    TEMPLATE_CLASSES  ${template_classes}
    SOURCES           ${sources}
    NOWRAP_CLASSES    ${nowrap_classes}
    HEADERS           ${headers}
    PRIVATE_TEMPLATES ${private_templates}
    PRIVATE_HEADERS   ${private_headers})

  vtk_module_link(VTK::ICP
    PRIVATE
    VTK::CommonDataModel)
endif()
  
