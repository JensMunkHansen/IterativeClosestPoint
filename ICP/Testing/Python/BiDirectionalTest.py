import os, sys

if os.name == 'posix':
  sys.path.insert(0, 3*"../" + "build/lib/python3.11/site-packages")
else:
  # On windows, we use the installed one
  sys.path.insert(0, 3*"../" + "install/bin/Lib/site-packages")

from icpmodules.util.scene_utils import (
    vtk_subfigs,
    vtk_show_points,
    vtk_show_gradients)
from icpmodules.util.misc import vtk_random_hills

# VTK includes

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2

from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkFiltersCore import  vtkPolyDataNormals
from vtkmodules.vtkRenderingCore import (
  vtkGlyph3DMapper,
  vtkActor,
  vtkPolyDataMapper,
  vtkRenderWindow,
  vtkRenderWindowInteractor,
  vtkRenderer)
from vtkmodules.vtkFiltersSources import vtkArrowSource
from vtkmodules.vtkCommonDataModel import (
    vtkPolyData,
    vtkDataObject,
    vtkDataSetAttributes)

from BiDirectional import BiDirectional

hills = vtk_random_hills()
hills.Update()
surfaceNormals = vtkPolyDataNormals()
surfaceNormals.SetInputConnection(hills.GetOutputPort())
surfaceNormals.Update()
sourceFilter = surfaceNormals

source = sourceFilter.GetOutput()

hills1 = vtk_random_hills(nHills=25)
hills1.Update()
target = hills1.GetOutput()

targetMapper = vtkPolyDataMapper()
targetMapper.SetInputData(target)
targetActor = vtkActor()
targetActor.SetMapper(targetMapper)

sourceMapper = vtkPolyDataMapper()
sourceMapper.SetInputConnection(sourceFilter.GetOutputPort())
sourceActor = vtkActor()
sourceActor.SetMapper(sourceMapper)

renWin, renderers = vtk_subfigs(2,2, sharecamera=True)

# Registration here

regTransform = vtkTransform()

icp = BiDirectional()
icp.SetSource(source)
icp.SetTarget(target)
# Find intersections
icp.InternalUpdate()


# Display results

if 1:
  regTransformPoly = vtkTransformPolyDataFilter()
  regTransformPoly.SetTransform(regTransform)
  regTransformPoly.SetInputConnection(sourceFilter.GetOutputPort())
  regTransformPoly.Update()
  
  regSourceMapper = vtkPolyDataMapper()
  regSourceMapper.SetInputConnection(regTransformPoly.GetOutputPort())
  regSourceActor = vtkActor()
  regSourceActor.SetMapper(regSourceMapper)
  regSourceActor.GetProperty().SetRepresentationToWireframe()
  regSourceActor.GetProperty().SetColor(1,0,0)
  
  renderers[0].AddActor(sourceActor)
  renderers[1].AddActor(targetActor)
  renderers[1].AddActor(regSourceActor)
  
  # Show points and gradients
  vtk_show_points(renderers[0], source)
  vtk_show_gradients(renderers[1],
                     icp.IntSourceToTarget.landmarks, # Issue
                     icp.IntSourceToTarget.normals)

if 1:
  invRegTransform = vtkTransform()
  invRegTransformPoly = vtkTransformPolyDataFilter()
  invRegTransformPoly.SetTransform(invRegTransform)
  invRegTransformPoly.SetInputData(target)
  invRegTransformPoly.Update()
  
  regTargetMapper = vtkPolyDataMapper()
  regTargetMapper.SetInputConnection(invRegTransformPoly.GetOutputPort())
  regTargetActor = vtkActor()
  regTargetActor.SetMapper(regTargetMapper)
  regTargetActor.GetProperty().SetRepresentationToWireframe()
  regTargetActor.GetProperty().SetColor(1,0,0)
  
  renderers[2].AddActor(sourceActor)
  renderers[2].AddActor(regTargetActor)
  renderers[3].AddActor(targetActor)
  
  vtk_show_points(renderers[3], target)
  vtk_show_gradients(renderers[2],
                     icp.IntTargetToSource.landmarks,
                     icp.IntTargetToSource.normals)


iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

renWin.Render()
renderers[1].ResetCamera()
iren.Start()
