vtk_add_test_cxx(
  vtkICPCxxTests tests
  TestPolyDataCorrespondenceFilter.cxx, NO_VALID
  TestKabsch.cxx, NO_VALID
)
vtk_test_cxx_executable(vtkICPCxxTests tests
  RENDERING_FACTORY
)

