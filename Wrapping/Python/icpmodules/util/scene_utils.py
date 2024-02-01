import re
import math

from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
  vtkStringArray,
  vtkMath,
  vtkUnsignedCharArray,
  vtkPoints,
  VTK_RGB)
from vtkmodules.vtkCommonDataModel import (
  vtkDataSet,
  vtkVector3d,
  vtkPolyData)
from vtkmodules.vtkCommonMath import (
  vtkMatrix3x3,
  vtkMatrix4x4,
  vtkQuaterniond)
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersCore import (
  vtkSmoothPolyDataFilter,
  vtkPolyDataNormals)
from vtkmodules.vtkFiltersSources import (
  vtkSphereSource,
  vtkArrowSource)
from vtkmodules.vtkInteractionWidgets import (
  vtkSliderRepresentation2D,
  vtkSliderWidget)
from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor
from vtkmodules.vtkRenderingCore import (
  vtkGlyph3DMapper,
  vtkRenderWindowInteractor,
  vtkRenderWindow,
  vtkCamera,
  vtkRenderer,
  vtkPolyDataMapper,
  vtkActor)

from vtkmodules.util.colors import white

def vtk_mapper_scalar_to_vertex_colors(actor, arrayname):

  # Do some renderering to trigger coloring
  renderer = vtkRenderer()
  renderer.AddActor(actor)
  renderWindow = vtkRenderWindow()
  renderWindow.SetOffScreenRendering(1)
  renderWindow.AddRenderer(renderer)
  renderWindow.Render()

  # Replace the colors with the mapped colors
  inputData = actor.GetMapper().GetInput()
  lut = actor.GetMapper().GetLookupTable()
  scalars = inputData.GetPointData().GetArray(arrayname)
  inputData.GetPointData().RemoveArray('RGB')
  colors = lut.MapScalars(scalars, 0, 0, VTK_RGB)
  if (colors.GetReferenceCount() == 2):
    # Leak in VTK
    colors.SetReferenceCount(1)
  colors.SetName('RGB')
  inputData.GetPointData().AddArray(colors)
  del renderWindow
  return inputData

def vtk_axes_to_transform(source, target):
  """
  Generate homegenous transform transforming origin and positive
  orientation defined by source[normal, first, origin] into target[normal, first, origin]
  """
  import numpy as np
  normal0 = source[0]
  first0 = source[1]
  origin0 = source[2]
  normal1 = target[0]
  first1 = target[1]
  origin1 = target[2]

  vec = vtkVector3d() # Axis of rotation
  vtkMath.Cross(normal0, normal1, vec)
  costheta = vtkMath.Dot(normal1, normal0)
  sintheta = vtkMath.Norm(vec)
  theta = np.arctan2(sintheta, costheta)

  if sintheta != 0.0:
    vec[0] = vec[0]/sintheta
    vec[1] = vec[1]/sintheta
    vec[2] = vec[2]/sintheta

  # Convert to Quaternion
  costheta = np.cos(0.5*theta)
  sintheta = np.sin(0.5*theta)
  quat0 = vtkQuaterniond(costheta, vec[0]*sintheta, vec[1]*sintheta, vec[2]*sintheta)

  newFirst = vtkVector3d()

  rot0 = np.ones((3,3),dtype=np.float64)
  vtkMath.QuaternionToMatrix3x3(quat0, rot0)

  if 1:
    # Can be performed using quaternions
    vtkMath.Multiply3x3(rot0, first0, newFirst)
  else:
    # Quaternion equivalent of the above line
    quatAxis0 = vtkQuaterniond(0.0, first0[0],
                                   first0[1],
                                   first0[2])
    quatAxisTmp = vtkQuaterniond()
    quatAxis1 = vtkQuaterniond()
    vtkMath.MultiplyQuaternion(quat0, quatAxis0, quatAxisTmp)
    vtkMath.MultiplyQuaternion(quatAxisTmp, quat0.Inverse(), quatAxis1)
    newFirst[0] = quatAxis1[1]
    newFirst[1] = quatAxis1[2]
    newFirst[2] = quatAxis1[3]

  # Rotate newFirst into first1
  vec = vtkVector3d() # Axis of rotation
  vtkMath.Cross(newFirst, first1, vec)
  costheta = vtkMath.Dot(first1, newFirst)
  sintheta = vtkMath.Norm(vec)
  theta = np.arctan2(sintheta, costheta)
  if sintheta != 0.0:
    vec[0] = vec[0]/sintheta
    vec[1] = vec[1]/sintheta
    vec[2] = vec[2]/sintheta

  # Convert to Quaternion
  costheta = np.cos(0.5*theta)
  sintheta = np.sin(0.5*theta)
  quat1 = vtkQuaterniond(costheta, vec[0]*sintheta, vec[1]*sintheta, vec[2]*sintheta)
  if 0:
    rot1 = np.ones((3,3),dtype=np.float64)
    vtkMath.QuaternionToMatrix3x3(quat1, rot1)
    rot = np.dot(rot1, rot0)
  else:
    # Quaternion equivalent of the above
    rot = np.ones((3,3),dtype=np.float64)
    quat2 = vtkQuaterniond()
    vtkMath.MultiplyQuaternion(quat1, quat0, quat2)
    vtkMath.QuaternionToMatrix3x3(quat2, rot)

  # Rotation
  mat = np.zeros((4,4), dtype=np.float64)
  mat[:3,:3] = rot
  mat[3,3] = 1.0

  # Translation
  tmp = vtkVector3d()
  vtkMath.Multiply3x3(rot, origin0, tmp)
  mat[:3,3] = np.array(origin1) - np.array(tmp)

  # Construct 4x4 matrix
  trans = vtkMatrix4x4()
  trans.DeepCopy(mat.flatten().tolist())

  return trans

def vtk_next_named_color(step=2):
  colors = vtkNamedColors()
  colorNames = colors.GetColorNames().split("\n")
  if not (vtk_next_named_color.iColor < len(colorNames)):
    vtk_next_named_color.iColor = 0
  colorName = colorNames[vtk_next_named_color.iColor]
  vtk_next_named_color.iColor = vtk_next_named_color.iColor + step
  return colors.GetColor3d(colorName)
vtk_next_named_color.iColor = 0

def vtk_origin_to_camera_matrix(cam):
  y = cam.GetViewUp()
  vtkMath.Normalize(y)
  
  position = cam.GetPosition()
  focalPoint = cam.GetFocalPoint()
  z = vtkVector3d()
  vtkMath.Subtract(focalPoint, position, z)
  vtkMath.Normalize(z)
  x = vtkVector3d()
  vtkMath.Cross(y, z, x)

  # Define the position
  position = cam.GetPosition()
  
  # Set the matrix of the vtkTransform using the vectors and position
  matrix = vtkMatrix4x4()
  matrix.SetElement(0, 0, x[0])
  matrix.SetElement(0, 1, y[0])
  matrix.SetElement(0, 2, z[0])
  matrix.SetElement(0, 3, position[0])
  matrix.SetElement(1, 0, x[1])
  matrix.SetElement(1, 1, y[1])
  matrix.SetElement(1, 2, z[1])
  matrix.SetElement(1, 3, position[1])
  matrix.SetElement(2, 0, x[2])
  matrix.SetElement(2, 1, y[2])
  matrix.SetElement(2, 2, z[2])
  matrix.SetElement(2, 3, position[2])
  matrix.SetElement(3, 3, 1)

  # ModelViewProjectionMatrix (MVP) 
  # ViewTransformMatrix (V) worldToCamera
  # ProjectionTransformMatrix (P) worldToScreen
  # M = VxP (actually it is PxV)

  # Create the transform - from origin to camera
  transform = vtkTransform()
  transform.SetMatrix(matrix)

  # Get the transform matrix (flipped z-axis and from camera to origin)
  view_matrix = cam.GetViewTransformMatrix()

  return transform

def vtk_create_axes_actor(polydata, color=white):
  cubeAxesActor = vtkCubeAxesActor()
  cubeAxesActor.SetUseTextActor3D(1)
  bounds = polydata.GetBounds()
  cubeAxesActor.SetBounds(bounds)
  cubeAxesActor.XAxisMinorTickVisibilityOff()
  cubeAxesActor.YAxisMinorTickVisibilityOff()
  cubeAxesActor.ZAxisMinorTickVisibilityOff()
  cubeAxesActor.SetFlyModeToStaticEdges()
  for i in range(3):
    cubeAxesActor.GetLabelTextProperty(i).SetColor(white)
    cubeAxesActor.GetTitleTextProperty(i).SetColor(white)
  cubeAxesActor.GetXAxesLinesProperty().SetColor(white)
  cubeAxesActor.GetYAxesLinesProperty().SetColor(white)
  cubeAxesActor.GetZAxesLinesProperty().SetColor(white)
  cubeAxesActor.GetProperty().SetColor(white)
  return cubeAxesActor

def vtk_wrap_surface(inputPolyData, normals=True):
  actor = vtkActor()
  bounds = inputPolyData.GetBounds()
  sphereSource = vtkSphereSource()
  
  radius = 1.1 * math.sqrt(sum([(i-j)*(i-j) for i,j in zip(bounds[1::2], bounds[0::2])]))
  
  sphereSource.SetRadius(radius)
  sphereSource.SetPhiResolution(80)
  sphereSource.SetThetaResolution(80)
  sphereSource.SetCenter(0.5*(bounds[0] + bounds[1]),
                         0.5*(bounds[2] + bounds[3]),
                         0.5*(bounds[4] + bounds[5]))

  smoothFilter = vtkSmoothPolyDataFilter()
  smoothFilter.SetInputConnection(0, sphereSource.GetOutputPort())
  smoothFilter.SetInputData(1, inputPolyData)

  normals = vtkPolyDataNormals()
  normals.SetInputConnection(smoothFilter.GetOutputPort())
  normals.SetFeatureAngle(60)
  normals.SplittingOff()
  normals.ComputeCellNormalsOn()
  normals.Update()

  mapper = vtkPolyDataMapper()
  mapper.SetInputConnection(normals.GetOutputPort())

  actor.SetMapper(mapper)

  return actor, normals.GetOutput()
  

def vtk_slider_2d(renderer, iren, _min=0, _max=10.0, title="Index"):
  SliderRepres = vtkSliderRepresentation2D()
  SliderRepres.SetRenderer(renderer)
  SliderRepres.SetMinimumValue(_min)
  SliderRepres.SetMaximumValue(_max)
  SliderRepres.SetValue(int((_min + _max) / 2))
  SliderRepres.SetTitleText(title)
  SliderRepres.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
  SliderRepres.GetPoint1Coordinate().SetValue(0.1, 0.05)
  SliderRepres.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
  SliderRepres.GetPoint2Coordinate().SetValue(0.4, 0.05)
  SliderRepres.SetSliderLength(0.02)
  SliderRepres.SetSliderWidth(0.03)
  SliderRepres.SetEndCapLength(0.01)
  SliderRepres.SetEndCapWidth(0.03)
  SliderRepres.SetTubeWidth(0.005)
  SliderRepres.SetLabelFormat("%3.2lf")
  SliderRepres.SetTitleHeight(0.02)
  SliderRepres.SetLabelHeight(0.02)

  # Slider widget
  SliderWidget = vtkSliderWidget()
  SliderWidget.SetInteractor(iren)
  SliderWidget.SetRepresentation(SliderRepres)
  SliderWidget.KeyPressActivationOff()
  SliderWidget.SetAnimationModeToAnimate()
  return SliderWidget

def vtk_show_points(renderer, points, color=[1,0,0]):
  if type(points) == vtkPolyData:
      actor = vtkActor()
      sphere = vtkSphereSource()
      sphere.SetRadius(0.08)
      pointMapper = vtkGlyph3DMapper()
      pointMapper.SetInputData(points)
      pointMapper.SetSourceConnection(sphere.GetOutputPort())
      pointMapper.ScalingOff()
      pointMapper.ScalarVisibilityOn()
      actor.SetMapper(pointMapper)
      actor.GetProperty().SetDiffuseColor(*color)
      renderer.AddActor(actor)
      return actor
  return None

def vtk_show_gradients(renderer, points, normals, color=(0,0,1)):
  if type(points) == vtkPoints:
    arrow_source = vtkArrowSource()
    # Create new polydata
    poly = vtkPolyData()
    poly.SetPoints(points)
    normals.SetName("GRADIENTS")
    poly.GetPointData().AddArray(normals)
    poly.GetPointData().SetActiveVectors("GRADIENTS")
    mapper = vtkGlyph3DMapper()
    mapper.SetScalarVisibility(0)
    mapper.SetInputData(poly)
    mapper.SetSourceConnection(arrow_source.GetOutputPort())
  
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(*color)
    renderer.AddActor(actor)
    return actor
  return None

def vtk_multiple_renderers(renderWindow, background=None, nrows=1, ncols=1, sharecamera=False):
  renderers = []
  camera = None
  if sharecamera:
    camera = vtkCamera()
  for irow in range(nrows-1,-1,-1):
    for icol in range(ncols):
      renderer = vtkRenderer()
      renderer.SetViewport(0.0 + icol*1.0/ncols,0.0 + irow*1.0/nrows,
                           (icol+1)*1.0/ncols,(irow+1)*1.0/nrows)
      if background is not None:
        renderer.SetBackground(background)
      if camera is not None:
        renderer.SetActiveCamera(camera)
      renderWindow.AddRenderer(renderer)
      renderers.append(renderer)
  return renderers

def vtk_subfigs(nrows=1, ncols=1, sharecamera=False):
  renderWindow = vtkRenderWindow()
  renderers = []
  camera = None
  if sharecamera:
    camera = vtkCamera()
  for irow in range(nrows-1,-1,-1):
    for icol in range(ncols):
      renderer = vtkRenderer()
      renderer.SetViewport(0.0 + icol*1.0/ncols,0.0 + irow*1.0/nrows,
                           (icol+1)*1.0/ncols,(irow+1)*1.0/nrows)
      if camera is not None:
        renderer.SetActiveCamera(camera)
      renderWindow.AddRenderer(renderer)
      renderers.append(renderer)
  return renderWindow, renderers

def vtk_acq_camera_info(camera:vtkCamera):
  stringBuffer = ""
  vector = camera.GetFocalPoint()
  stringBuffer +="Camera:FocalPoint " + str(vector[0]) + ", " + str(vector[1]) + ", " + str(vector[2]) + "\n"
  vector = camera.GetPosition()
  stringBuffer +="Camera:Position " + str(vector[0]) + ", " + str(vector[1]) + ", " + str(vector[2]) + "\n"
  vector = camera.GetViewUp()
  stringBuffer +="Camera:ViewUp " + str(vector[0]) + ", " + str(vector[1]) + ", " + str(vector[2]) + "\n"
  scalar = camera.GetViewAngle()
  stringBuffer +="Camera:ViewAngle " + str(scalar) + "\n";
  vector = camera.GetClippingRange()
  stringBuffer +="Camera:ClippingRange " + str(vector[0]) + ", " + str(vector[1])
  return stringBuffer

def vtk_save_scene_to_field_data(data:vtkDataSet, camera:vtkCamera):
  stringBuffer = vtk_acq_camera_info(camera)
  cameraArray = vtkStringArray()
  cameraArray.SetNumberOfValues(1)
  cameraArray.SetValue(0, stringBuffer)
  cameraArray.SetName("Camera")
  data.GetFieldData().AddArray(cameraArray)

def vtk_restore_scene_from_field_data(data:vtkDataSet, camera:vtkCamera):
  lines = None
  # Get the saved camera information from the field data
  if (vtkStringArray.SafeDownCast(
        data.GetFieldData().GetAbstractArray("Camera"))):
    lines = vtkStringArray.SafeDownCast(data.GetFieldData().GetAbstractArray("Camera")).GetValue(0)
  else:
    return

  reCP = "^Camera:Position (.+)"
  reCFP = "^Camera:FocalPoint (.+)"
  reCVU = "^Camera:ViewUp (.+)"
  reCVA = "^Camera:ViewAngle (.+)"
  reCCR = "^Camera:ClippingRange (.+)"
  floatNumber = "[^0-9\\.\\-]*([0-9e\\.\\-]*[^,])[^0-9\\.\\-]*([0-9e\\.\\-]*[^,])[^0-9\\.\\-]*([0-9e\\.\\-]*[^,])"
  floatScalar = "[^0-9\\.\\-]*([0-9\\.\\-e]*[^,])"

  for line in lines.split("\n"):
    if match := re.search(reCFP, line):
      rest = match.groups(0)
      if match := re.search(floatNumber, rest[0]):
        x = float(match.groups(0)[0])
        y = float(match.groups(0)[1])
        z = float(match.groups(0)[2])
        camera.SetFocalPoint(x,y,z)
    elif match := re.search(reCP, line):
      rest = match.groups(0)[0]
      if match := re.search(floatNumber, rest):
        x = float(match.groups(0)[0])
        y = float(match.groups(0)[1])
        z = float(match.groups(0)[2])
        camera.SetPosition(x, y, z)
    elif match := re.search(reCVU, line):
      rest = match.groups(0)[0]
      if match := re.search(floatNumber, rest):
        x = float(match.groups(0)[0])
        y = float(match.groups(0)[1])
        z = float(match.groups(0)[2])
        camera.SetViewUp(x, y ,z)
    elif match := re.search(reCVA, line):
      rest = match.groups(0)[0]
      if match := re.search(floatScalar, rest):
        camera.SetViewAngle(float(match.groups(0)[0]))
    elif match := re.search(reCCR, line):
      rest = match.groups(0)[0]
      if match := re.search(floatNumber, rest):
        camera.SetClippingRange(float(match.groups(0)[0]),
                                float(match.groups(0)[1]))
    else:
      print("Unrecognized line: " + line)

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
