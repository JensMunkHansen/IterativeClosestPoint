#pragma once

#include "vtkMultiMeshAlignment.h"
#include <list>

class vtkMultiMeshAlignment::vtkInternals
{
public:
  struct PolyDataPair
  {
    // Moving
    vtkPolyData* Source;
    // Reference
    vtkPolyData* Target;
  };
  vtkMultiMeshAlignment* Parent;
  vtkInternals(vtkMultiMeshAlignment* parent);
  ~vtkInternals();
};
