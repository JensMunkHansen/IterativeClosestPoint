import sys
import os
import socket

# Tuesday
# Jacobian
#
# Loop with line-search
from icpmodules.vtkICP import vtkImplicitPolyDataIntersection

from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline

from vtkmodules.vtkCommonCore import (
  vtkDoubleArray,
  vtkFloatArray,
  vtkMath,
  vtkPoints,
  reference,
)

from vtkmodules.vtkCommonDataModel import (
  vtkDataObject,
  vtkDataSet,
  vtkFieldData,
  vtkPolyData,
  vtkCellLocator,
  vtkGenericCell,
)

from vtkmodules.vtkFiltersCore import vtkTriangleFilter
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkCommonTransforms import vtkTransform

from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtk.numpy_interface import dataset_adapter as dsa

import numpy as np

from euler import euler2rot, euler2rotd

class Intersections():
  pass

# From source to target.
class BiDirectional(object):
  def __init__(self):
    self._source = None
    self._target = None
    self._nMaxLandmarks = 2000
    self.VnLine = np.r_[0.0, 0.0, 1.0]
    self.maxDist = 16.3
    self.Offset = np.r_[0.0, 0.0, 0.0]
    self.Angle = np.r_[0.0, 0.0, 0.0]
    self.CreateDefaultLocators()
    self.SourceSubscanToWorld = vtkTransform()
    self.TargetSubscanToWorld = vtkTransform()
    pass
  def CreateDefaultLocators(self):
    self._locatorTarget = vtkImplicitPolyDataIntersection()
    self._locatorSource = vtkImplicitPolyDataIntersection()
  def PrepareData(self, data):
    # Keep only polygons
    tri = vtkTriangleFilter()
    tri.PassVertsOff()
    tri.PassLinesOff()
    tri.SetInputData(data)
    tri.Update()
    return tri.GetOutput()
  def SetSource(self, source):
    if (source != self._source):
      self._source = self.PrepareData(source)
      self._locatorSource.SetInput(self._source)
  def SetTarget(self, target):
    if (target != self._target):
      self._target = self.PrepareData(target)
      self._locatorTarget.SetInput(self._target)      
  def FindIntersections(self, direction, transform):
    """
    Transform used for transforming point and direction of source(moving)
    """
    halfVnLine = 0.5 * self.maxDist * self.VnLine

    if direction == 0:
      source = self._source
      locator = self._locatorTarget
    else:
      source = self._target
      locator = self._locatorSource

    nSourcePoints = source.GetNumberOfPoints()
    step = 1
    if (step > self._nMaxLandmarks):
      step = int(nSourcePoints / self._nMaxLandmarks)
    nSourcePoints = int(nSourcePoints / step)
    
    signedDistance = reference(0.0)
    hitPoint = [0,0,0]
    gradient = [0,0,0]
    lineP0 = [0,0,0]
    lineP1 = [0,0,0]

    points = vtkPoints()
    landmarks = vtkPoints()
    points.Initialize()
    landmarks.Initialize()
    transform.InternalTransformVector(halfVnLine, halfVnLine)
    vnLine = np.r_[0.0, 0.0, 0.0]
    transform.InternalTransformVector(self.VnLine, vnLine)
    normals = vtkDoubleArray()
    normals.SetNumberOfComponents(3)
    weightsum = 0.0
    j = 0
    for i in range(nSourcePoints):
      # Line in world coordinates
      sourcePoint = source.GetPoint(j)
      transform.InternalTransformPoint(sourcePoint, sourcePoint)
      vtkMath.Subtract(sourcePoint, halfVnLine, lineP0)
      vtkMath.Add(sourcePoint, halfVnLine, lineP1)

      # Compute intersection in world coordinates
      nIntersections = locator.IntersectWithLine(lineP0, lineP1, signedDistance, hitPoint, gradient)
      if (nIntersections > 0):
        # Same coordinates system: 3S (refSubscan), JMH (world)
        landmarks.InsertNextPoint(hitPoint)
        points.InsertNextPoint(sourcePoint) 
        normals.InsertNextTuple(gradient)
        weightsum = weightsum + 1.0
      j = j + step
    intersections = Intersections()
    intersections.points = points
    intersections.landmarks = landmarks
    intersections.normals = normals
    intersections.weightsum = weightsum
    intersections.weights = np.ones(points.GetNumberOfPoints())
    intersections.vnLine = vnLine
    return intersections
  def InternalUpdate(self, transform:vtkTransform = None):
    if (self._target is None or self._target.GetNumberOfPoints() < 1):
      return
    if (self._source is None or self._source.GetNumberOfPoints() < 1):
      return

    if transform is None:
      accumulatedTransform = vtkTransform()
      accumulatedTransform.Identity()
    # From source to world -> target (world)
    accumulatedTransform.PostMultiply()
    #accumulatedTransform.Concatenate(self.SourceSubscanToWorld)
    accumulatedTransform.Update()

    # Source -> Target
    self.IntSourceToTarget = self.FindIntersections(0, accumulatedTransform)
    # Target -> Source
    self.IntTargetToSource = self.FindIntersections(1, accumulatedTransform)

    self.Evaluate()

  def Evaluate(self):
    dX, dY, dZ = euler2rotd(self.Angle, conv='xyz', intrinsic=False)
    rot = euler2rot(self.Angle, conv='xyz', intrinsice=False)
    costNorm = 0.0
    costInv = 0.0
    costNorm, jTjNorm, jTrNorm = self.EvaluateJacobian(self.IntSourceToTarget, self.Offset, rot, dX, dY, dZ)
    costInv, jTjInv, jTrInv = self.EvaluateJacobianInv(self.IntSourceToTarget, self.Offset, rot, dX, dY, dZ)
    #jTj = jTjNorm + jTjInv
    #jTr = jTrNorm + jTrInv
    return costNorm# + costInv, jTj, jTr
    
  def EvaluateJacobian(self, intersections, translation, rot, dX, dY, dZ):
    points = vtk_to_numpy(intersections.points.GetData())
    landmarks = vtk_to_numpy(intersections.landmarks.GetData())
    normals = vtk_to_numpy(intersections.normals)
    vnLine = intersections.vnLine

    result = 0.0
    jTj = np.zeros((6,6))
    jTr = np.zeros(6)
    maxSqrDist = self.maxDist * self.maxDist
    vnRay = np.dot(rot, vnLine)
    print(vnRay)
    for i in range(len(points)):
      ptLine = landmarks[i]
      ptRay = np.dot(rot,ptLine) + translation
      # Line-plane distance
      normal = normals[i]
      denominator = 1.0 / np.dot(normal, vnRay)
      numerator = np.dot((landmarks[i] - ptRay), normal)
      l = numerator * denominator
      l2 = l * l
      if (l2 < maxSqrDist):
        dldTranslation = -denominator * normal
        factor = (-l * vnLine - landmarks[i]) * denominator
        dldR = np.r_[np.dot(np.dot(dX, factor), normal),
                     np.dot(np.dot(dY, factor), normal),
                     np.dot(np.dot(dZ, factor), normal)]
        jac = np.r_[dldTranslation, dldR]
        jTj = jTj + intersections.weights[i] * np.outer(jac, jac)
        jTr = jTr + jac * intersections.weights[i]
        result = result + l2 * intersections.weights[i]
      else:
        result = result + maxSqrDist * intersections.weights[i]
    return result, jTj, jTr
  def EvaluateJacobianInv(self, intersections, translation, rot, dX, dY, dZ):
    points = vtk_to_numpy(intersections.points.GetData())
    landmarks = vtk_to_numpy(intersections.landmarks.GetData())
    normals = vtk_to_numpy(intersections.normals)
    vnLine = intersections.vnLine

    result = 0.0
    jTj = np.zeros((6,6))
    jTr = np.zeros(6)
    maxSqrDist = self.maxDist * self.maxDist

    for i in range(len(points)):
      ptPlaneTransformed = np.dot(rot, landmarks[i]) + translation
      vnPlaneTransformed = np.dot(rot, normals[i])
      ptDiff = ptPlaneTransformed - landmarks[i]
      denominator = 1.0 / np.dot(vnLine, vnPlaneTransformed)
      l = denominator * np.dot(ptDiff, vnPlaneTransformed)
      l2 = l*l
      if l2 < maxSqrDist:
          # dLdTranslation is the derivative of l with respect to translation
          dLdTranslation = denominator * vnPlaneTransformed;
          factor = (ptDiff - l * vnLine) * denominator;
          # dLdR is the derivative of l with respect to the euler angles through the chain rule,
          dLdR = np.r_[
            np.dot(dLdTranslation, np.dot(dX, landmarks[i])) + np.dot(factor, np.dot(dX, normals[i])),
            np.dot(dLdTranslation, np.dot(dY, landmarks[i])) + np.dot(factor, np.dot(dY, normals[i])),
            np.dot(dLdTranslation, np.dot(dZ, landmarks[i])) + np.dot(factor, np.dot(dZ, normals[i]))]
          jac = np.r_[dLdTranslation, dLdR]
          jTj = jTj + np.outer(jac,jac) * intersections.weights[i]
          jTr = jac + l * intersections.weights[i]
          result = result + l2 * intersections.weights[i]
      else:
          result = result + maxSqrDist * intersections.weights[i];
    return result, jTj, jTr
  def ComputeIndices(self, source):
    nPoints = source.GetNumberOfPoints()
    if (nPoints > self._nMaxLandmarks):
      step = int(nPoints / self._nMaxLandmarks)
    else:
      step = 1
    return np.r_[0:nPoints:step]
# points, direction: world-to-subscan
# line: point + direction
# Get vertices around intersection
# Get normal of facet
# Normal, point
# Weights (sum)
# Line direction in subscan space
    
  def BuildLandmarks(self, source):
    if (source.GetNumberOfPoints() > self._nMaxLandmarks):
      step = int(source.GetNumberOfPoints() / self._nMaxLandmarks)
    else:
      step = 1
    nSourcePoints = int(source.GetNumberOfPoints() / step)

    if (source == self._source):
      locator = self._locatorTarget
    else:
      locator = self._locatorSource
    j = 0
    signedDistance = reference(0.0)
    hitPoint = [0,0,0]
    gradient = [0,0,0]

    points = vtkPoints()
    landmarks = vtkPoints()
    for i in range(nSourcePoints):
      lineP0 = [0,0,0]
      vtkMath.Subtract(source.GetPoint(j), (0,0,0.5), lineP0)
      lineP1 = [0,0,0]
      vtkMath.Add(source.GetPoint(j), (0,0,0.5), lineP1)
      nIntersections = locator.IntersectWithLine(lineP0, lineP1, signedDistance, hitPoint, gradient)
      if (nIntersections > 0):
        landmarks.InsertNextPoint(hitPoint)
        points.InsertNextPoint(source.GetPoint(j))
      j = j + step
    return points, landmarks
    
  def UpdateLandmarks(self):
    pass
  
# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
    
