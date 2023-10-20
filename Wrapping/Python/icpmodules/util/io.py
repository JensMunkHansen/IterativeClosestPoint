#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from vtkmodules.vtkIOLegacy import vtkPolyDataReader
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.vtkIOGeometry import (
  vtkBYUReader,
  vtkSTLReader,
  vtkOBJReader
)
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkFiltersCore import (
  vtkAppendPolyData
)
class vtkPolyDataReaderFactory:
  def __new__(cls):
    raise TypeError("Static classes cannot be instantiated")
  @staticmethod
  def CreateReader(inputFileName):
    fileName, fileExtension = os.path.splitext(inputFileName)
    if fileExtension in ['.ply', '.PLY']:
      reader = vtkPLYReader()
    elif fileExtension in ['.vtp', '.VTP']:
      reader = vtkXMLPolyDataReader()
    elif fileExtension in ['.stl', '.STL']:
      reader = vtkSTLReader()
    elif fileExtension in ['.vtk', '.VTK']:
      reader = vtkPolyDataReader()
    elif fileExtension in ['.g', '.G']:
      reader = vtkBYUReader()
    elif fileExtension in ['.obj', '.OBJ']:
      reader = vtkOBJReader()
    else:
      raise NameError("Filetype not supported")
    reader.SetFileName(inputFileName)
    return reader
  @staticmethod
  def SupportedExtensions():
    return [".ply", ".vtp", ".stl", ".vtk", ".g", ".obj"]
  @staticmethod
  def CreateReaders(inputFileNames, append=False):
    readers = []
    for inputFileName in inputFileNames:
      readers.append(vtkPolyDataReaderFactory.CreateReader(inputFileName))
    if not append:
      return readers
    else:
      append = vtkAppendPolyData()
      for reader in readers:
        append.AddInputConnection(reader.GetOutputPort())
      return append

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
