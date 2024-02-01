from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkFiltersSources import vtkParametricFunctionSource  
from vtkmodules.vtkCommonComputationalGeometry import vtkParametricRandomHills

def vtk_random_hills(nHills = 20):
  sourceFcn = vtkParametricRandomHills()
  sourceFcn.AllowRandomGenerationOn()
  sourceFcn.SetRandomSeed(1)
  sourceFcn.SetNumberOfHills(nHills)
  source = vtkParametricFunctionSource()
  source.SetParametricFunction(sourceFcn)
  source.SetUResolution(51)
  source.SetVResolution(51)
  source.SetWResolution(51)
  source.Update()
  return source
