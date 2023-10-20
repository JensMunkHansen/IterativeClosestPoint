set(sources
  vtkPolyDataCorrespondenceFilter.cxx
  vtkICP.cxx
  vtkImplicitPolyDataDistance2.cxx
)

set(private_headers
  icp/icp.hpp
  icp/point2point.hpp
  icp/point2plane.hpp
)
set(nowrap_classes
  icp/icp
  icp/point2plane
  icp/point2point)

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
  
