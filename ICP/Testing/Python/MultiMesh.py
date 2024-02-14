import sys
sys.path.insert(0, 3*"../" + "build/lib/python3.11/site-packages")
import icp

meshAlignment = icp.vtkMultiMeshAlignment()

from vtkmodules.vtkCommonDataModel import vtkPolyDataCollection
from vtkmodules.vtkCommonTransforms import vtkTransformCollection, vtkTransform
from vtkmodules.vtkFiltersSources import vtkSphereSource

collection = vtkPolyDataCollection()
transforms = vtkTransformCollection()
for i in range(4):
    sphereSource = vtkSphereSource()
    sphereSource.Update()
    collection.AddItem(sphereSource.GetOutput())
    transform = vtkTransform()
    transform.Translate(0.1, 0.0, 0.0)
    transform.Update()
    transforms.AddItem(transform)
meshAlignment.SetPolyDataCollection(collection)
meshAlignment.SetTransformCollection(transforms)
allGood = True
for i in range(4):
    for j in range(4):
        if j != i:
            allGood = allGood and meshAlignment.PairPolyData(collection.GetItemAsObject(i),
                                                             collection.GetItemAsObject(j), transforms.GetItemAsObject(0))
