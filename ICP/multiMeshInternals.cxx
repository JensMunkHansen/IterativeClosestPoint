#include "multiMeshInternals.h"
#include "vtkMultiMeshAlignment.h"

vtkMultiMeshAlignment::vtkInternals::vtkInternals(vtkMultiMeshAlignment* parent)
  : Parent(parent)
{
}

vtkMultiMeshAlignment::vtkInternals::~vtkInternals() {}
