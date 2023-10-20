set(LIB_SRCS
  vtkPolyDataCorrespondenceFilter.cxx
  vtkICP.cxx
  vtkImplicitPolyDataDistance2.cxx
)

vtk_module_add_module(VTK::ICP
  SOURCES ${LIB_SRCS}
  HEADERS ${LIB_SRC_HDRS} ${LIB_HDRS})  
