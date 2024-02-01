#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
# Improved version
filedir = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0, os.path.join(filedir, '../../../build/lib/python3.11/site-packages'))

from icpmodules.vtkICP import vtkICP
from icpmodules.util.io import vtkPolyDataReaderFactory

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
    VTK_DOUBLE_MAX,
    vtkPoints
)
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkCommonDataModel import (
    vtkIterativeClosestPointTransform,
    vtkPolyData
)
from vtkmodules.vtkCommonTransforms import (
    vtkLandmarkTransform,
    vtkTransform
)
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkFiltersCore import vtkPolyDataNormals

from vtkmodules.vtkFiltersModeling import vtkHausdorffDistancePointSetFilter
from vtkmodules.vtkIOGeometry import (
    vtkBYUReader,
    vtkOBJReader,
    vtkSTLReader
)
from vtkmodules.vtkIOLegacy import (
    vtkPolyDataReader,
    vtkPolyDataWriter
    )
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
from vtkmodules.vtkCommonCore import vtkSMPTools

smp = vtkSMPTools()
smp.SetBackend("TBB")

from timeit import default_timer as timer

def get_program_parameters():
    import argparse
    description = 'Align a vtkPolyData to a transformed one, e.g. Bunny.vtp.'
    epilogue = '''
    The Bunny.vtp can be found in the data directory.
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('src_fn', help='The polydata source file name.')

    args = parser.parse_args()

    return args.src_fn

colors = vtkNamedColors()

src_fn = get_program_parameters()
print('Loading source:', src_fn)
reader = vtkPolyDataReaderFactory.CreateReader(src_fn)
reader.Update()
target_polydata = reader.GetOutput()

# Compute normals for the target
normalsFilter = vtkPolyDataNormals()
normalsFilter.SetInputData(target_polydata);
normalsFilter.ComputeCellNormalsOn()
normalsFilter.Update()

target_polydata2 = normalsFilter.GetOutput()

original_source_polydata = vtkPolyData()
original_source_polydata.DeepCopy(target_polydata)

# Transform the source (subject for registration)
trnf = vtkTransform()
trnf.Translate(0.0,0.0, 0.3)
trnf.RotateY(15)

tpd = vtkTransformPolyDataFilter()
tpd.SetTransform(trnf)
tpd.SetInputData(original_source_polydata)
tpd.Update()
source_polydata = tpd.GetOutput()

# Setup render windows
renderer = vtkRenderer()
render_window = vtkRenderWindow()
render_window.AddRenderer(renderer)
interactor = vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

def ComputeHausdorffDistance(target_polydata, source_polydata):
    distance = vtkHausdorffDistancePointSetFilter()
    distance.SetInputData(0, target_polydata)
    distance.SetInputData(1, source_polydata)
    distance.Update()
    hausdorff = distance.GetOutput(0).GetFieldData().GetArray('HausdorffDistance').GetComponent(0, 0)
    return hausdorff

distance_before_align = ComputeHausdorffDistance(target_polydata, source_polydata)
print("Hausdorff distance before alignment: %f" % (distance_before_align))

# Original source polydata
original_source_mapper = vtkDataSetMapper()
original_source_mapper.SetInputData(original_source_polydata) # White
original_source_mapper.ScalarVisibilityOff()

original_source_actor = vtkActor()
original_source_actor.SetMapper(original_source_mapper)
original_source_actor.GetProperty().SetOpacity(0.6)
original_source_actor.GetProperty().SetDiffuseColor(
    colors.GetColor3d('White'))
renderer.AddActor(original_source_actor)

# Target polydata
target_mapper = vtkDataSetMapper()
target_mapper.SetInputData(tpd.GetOutput())
target_mapper.ScalarVisibilityOff()

target_actor = vtkActor()
target_actor.SetMapper(target_mapper)
target_actor.GetProperty().SetDiffuseColor(
    colors.GetColor3d('Tomato'))
renderer.AddActor(target_actor)

running_source_polydata = vtkPolyData()
running_source_polydata.DeepCopy(source_polydata)

sourceTransform = vtkTransform()
sourceTransform.Identity()
sourceTransform.PostMultiply()

tpd1 = vtkTransformPolyDataFilter()
tpd1.SetInputData(source_polydata)
tpd1.SetTransform(sourceTransform)
tpd1.Update()

# Source (registered to target)
regActor = vtkActor()
regMapper = vtkDataSetMapper()
regMapper.SetInputData(tpd1.GetOutput())
regMapper.ScalarVisibilityOff()

regActor = vtkActor()
regActor.SetMapper(regMapper)
regActor.GetProperty().SetOpacity(0.6)
regActor.GetProperty().SetDiffuseColor(
    colors.GetColor3d('Blue'))
renderer.AddActor(regActor)

# Standard iterative closest point from VTK
icp = vtkIterativeClosestPointTransform()
icp.GetLandmarkTransform().SetModeToRigidBody()
icp.SetMaximumNumberOfIterations(1)
icp.StartByMatchingCentroidsOff()
icp.SetMaximumNumberOfLandmarks(100)
icp.SetMaximumMeanDistance(.00001)
icp.CheckMeanDistanceOn()
icp.SetTarget(target_polydata)
icp.SetSource(running_source_polydata)

tpd2 = vtkTransformPolyDataFilter()

# Registration using point-to-plane metric
icpr = vtkICP()

icpr.CheckMeanDistanceOn()
icpr.SetMaximumMeanDistance(0.001)
icpr.SetMetricToPointToPoint()
#icpr.SetMetricToPointToPlane() # BUG
icpr.SetMetricToPlaneToPlane()
icpr.SetMaximumNumberOfIterations(100)
icpr.SetMaximumNumberOfLandmarks(100)
icpr.SetSource(running_source_polydata)
icpr.SetTarget(target_polydata2)
start = timer()
icpr.Update()
elapsed = timer() - start
print("Elapsed seconds: %f" % (elapsed))
print("Number of iterations: %d" % (icpr.GetNumberOfIterations()))
print("Mean distance: %f" % (icpr.GetMeanDistance()))

# Animation registration using standard iterative closest point from VTK
global count
count = 0;
def UpdateVTKRegistration(interactor, event):
  global count
  icp.Update()
  # Incremental transform of source 
  mat = vtkMatrix4x4()
  mat.DeepCopy(icp.GetMatrix())
  sourceTransform.Concatenate(mat)

  # Update display
  start = timer()
  tpd1.Update()
  elapsed = timer() - start
  regMapper.Update()
  render_window.Render()
  tpd2.SetTransform(icp)
  tpd2.SetInputData(running_source_polydata)
  tpd2.Update()
    
  # Replace source
  running_source_polydata.DeepCopy(tpd2.GetOutput())
  icp.SetSource(running_source_polydata)
  icp.Modified()

  # Check hausdorf distance
  distance_after_align = ComputeHausdorffDistance(target_polydata, running_source_polydata)
  render_window.Render()
  if (count > 1000 or distance_after_align < 1e-3):
    interactor.DestroyTimer(timerId)
    print("Number of iteration: %d" % (count))
    print("Hausdorff distance: %f" % (distance_after_align))
    return
  count = count + 1

interactor.Initialize()
interactor.AddObserver("TimerEvent", UpdateVTKRegistration)
timerId = interactor.CreateRepeatingTimer(100)

render_window.Render()

interactor.Start()

