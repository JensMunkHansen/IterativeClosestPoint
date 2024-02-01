/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImplicitPolyDataDistance.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImplicitPolyDataIntersection.h"

#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkSmartPointer.h"
#include "vtkTriangleFilter.h"

vtkStandardNewMacro(vtkImplicitPolyDataIntersection);

//------------------------------------------------------------------------------
vtkImplicitPolyDataIntersection::vtkImplicitPolyDataIntersection()
{
  this->NoClosestPoint[0] = 0.0;
  this->NoClosestPoint[1] = 0.0;
  this->NoClosestPoint[2] = 0.0;

  this->NoGradient[0] = 0.0;
  this->NoGradient[1] = 0.0;
  this->NoGradient[2] = 1.0;

  this->NoValue = 0.0;

  this->Input = nullptr;
  this->Locator = nullptr;
  this->Tolerance = 1e-12;
}

//------------------------------------------------------------------------------
void vtkImplicitPolyDataIntersection::SetInput(vtkPolyData* input)
{
  if (this->Input != input)
  {
    // Use a vtkTriangleFilter on the polydata input.
    // This is done to filter out lines and vertices to leave only
    // polygons which are required by this algorithm for cell normals.
    vtkNew<vtkTriangleFilter> triangleFilter;
    triangleFilter->PassVertsOff();
    triangleFilter->PassLinesOff();

    triangleFilter->SetInputData(input);
    triangleFilter->Update();

    this->Input = triangleFilter->GetOutput();

    this->Input->BuildLinks();
    this->NoValue = this->Input->GetLength();

    this->CreateDefaultLocator();
    this->Locator->SetDataSet(this->Input);
    this->Locator->SetTolerance(this->Tolerance);
    this->Locator->SetNumberOfCellsPerBucket(10);
    this->Locator->CacheCellBoundsOn();
    this->Locator->AutomaticOn();
    this->Locator->BuildLocator();
  }
}

//------------------------------------------------------------------------------
vtkMTimeType vtkImplicitPolyDataIntersection::GetMTime()
{
  vtkMTimeType mTime = this->vtkObject::GetMTime();
  vtkMTimeType InputMTime;

  if (this->Input != nullptr)
  {
    InputMTime = this->Input->GetMTime();
    mTime = (InputMTime > mTime ? InputMTime : mTime);
  }

  return mTime;
}

//------------------------------------------------------------------------------
vtkImplicitPolyDataIntersection::~vtkImplicitPolyDataIntersection()
{
  if (this->Locator)
  {
    this->Locator->UnRegister(this);
    this->Locator = nullptr;
  }
}

//------------------------------------------------------------------------------
void vtkImplicitPolyDataIntersection::CreateDefaultLocator()
{
  if (this->Locator == nullptr)
  {
    this->Locator = vtkCellLocator::New();
  }
}

int vtkImplicitPolyDataIntersection::IntersectWithLine(const double p1[3], const double p2[3],
  double& ret, double closestPoint[3], double g[3], vtkDataObject::AttributeTypes* type)
{
  // Set defaults
  ret = this->NoValue;

  for (int i = 0; i < 3; i++)
  {
    g[i] = this->NoGradient[i];
  }

  for (int i = 0; i < 3; i++)
  {
    closestPoint[i] = this->NoClosestPoint[i];
  }

  // See if data set with polygons has been specified
  if (this->Input == nullptr || Input->GetNumberOfCells() == 0)
  {
    vtkErrorMacro(<< "No polygons to intersect!");
    return 0;
  }

  double p[3];
  vtkIdType cellId;
  int subId;

  vtkDataArray* cnorms = nullptr;
  if (this->Input->GetCellData() && this->Input->GetCellData()->GetNormals())
  {
    // We have cell normals
    cnorms = this->Input->GetCellData()->GetNormals();
  }
  double pCoords[3];
  double t = 0.0;

  auto cell = this->TLCell.Local();
  int intersectionFound =
    this->Locator->IntersectWithLine(p1, p2, this->Tolerance, t, p, pCoords, subId, cellId, cell);
  if (!intersectionFound)
  {
    return intersectionFound;
  }
  else
  {
    ret = std::sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) * t;
  }

  if (cellId != -1) // point located
  {
    int count = 0;

    // grad = (point - p1) / dist
    for (int i = 0; i < 3; i++)
    {
      g[i] = (p[i] - p1[i]) / (ret == 0. ? 1. : ret);
    }

    double dist2, weights[3], pcoords[3], awnorm[3] = { 0, 0, 0 };
    cell->EvaluatePosition(p, closestPoint, subId, pcoords, dist2, weights);

    auto idList = this->TLCellIds.Local();
    for (int i = 0; i < 3; i++)
    {
      count += (std::abs(weights[i]) < this->Tolerance ? 1 : 0);
    }
    // Face case - weights contains no 0s
    if (count == 0)
    {
      *type = vtkDataObject::CELL;
      // Compute face normal.
      if (cnorms)
      {
        // We have cell normals
        cnorms->GetTuple(cellId, awnorm);
      }
      else
      {
        // Compute normal for the given cell
        vtkPolygon::ComputeNormal(cell->Points, awnorm);
      }
    }
    // Edge case - weights contain one 0
    else if (count == 1)
    {
      *type = vtkDataObject::EDGE;

      // ... edge ... get two adjacent faces, compute average normal
      int a = -1, b = -1;
      for (int edge = 0; edge < 3; edge++)
      {
        if (std::abs(weights[edge]) < this->Tolerance)
        {
          a = cell->GetPointId((edge + 1) % 3);
          b = cell->GetPointId((edge + 2) % 3);
          break;
        }
      }

      if (a == -1)
      {
        vtkErrorMacro(<< "Could not find edge when closest point is "
                      << "expected to be on an edge.");
        return this->NoValue;
      }

      // The first argument is the cell ID. We pass a bogus cell ID so that
      // all face IDs attached to the edge are returned in the idList.
      this->Input->GetCellEdgeNeighbors(VTK_ID_MAX, a, b, idList);
      for (int i = 0; i < idList->GetNumberOfIds(); i++)
      {
        double norm[3];
        if (cnorms)
        {
          cnorms->GetTuple(idList->GetId(i), norm);
        }
        else
        {
          this->Input->GetCell(idList->GetId(i), cell);
          vtkPolygon::ComputeNormal(cell->GetPoints(), norm);
        }
        awnorm[0] += norm[0];
        awnorm[1] += norm[1];
        awnorm[2] += norm[2];
      }
      vtkMath::Normalize(awnorm);
    }

    // Vertex case - weights contain two 0s
    else if (count == 2)
    {
      *type = vtkDataObject::VERTEX;
      // ... vertex ... this is the expensive case, get all adjacent
      // faces and compute sum(a_i * n_i) Angle-Weighted Pseudo
      // Normals, J. Andreas Baerentzen and Henrik Aanaes
      int a = -1;
      for (int i = 0; i < 3; i++)
      {
        if (std::abs(weights[i]) > this->Tolerance)
        {
          //          a = cell->PointIds->GetId(i);
          a = cell->GetPointId(i);
        }
      }

      if (a == -1)
      {
        vtkErrorMacro(<< "Could not find point when closest point is "
                      << "expected to be a point.");
        return this->NoValue;
      }

      this->Input->GetPointCells(a, idList);
      for (int i = 0; i < idList->GetNumberOfIds(); i++)
      {
        double norm[3];
        this->Input->GetCell(idList->GetId(i), cell);
        if (cnorms)
        {
          cnorms->GetTuple(idList->GetId(i), norm);
        }
        else
        {
          vtkPolygon::ComputeNormal(cell->GetPoints(), norm);
        }

        // Compute angle at point a
        int b = cell->GetPointId(0);
        int c = cell->GetPointId(1);
        if (a == b)
        {
          b = cell->GetPointId(2);
        }
        else if (a == c)
        {
          c = cell->GetPointId(2);
        }
        double pa[3], pb[3], pc[3];
        this->Input->GetPoint(a, pa);
        this->Input->GetPoint(b, pb);
        this->Input->GetPoint(c, pc);
        for (int j = 0; j < 3; j++)
        {
          pb[j] -= pa[j];
          pc[j] -= pa[j];
        }
        vtkMath::Normalize(pb);
        vtkMath::Normalize(pc);
        double alpha = std::acos(vtkMath::Dot(pb, pc));
        awnorm[0] += alpha * norm[0];
        awnorm[1] += alpha * norm[1];
        awnorm[2] += alpha * norm[2];
      }
      vtkMath::Normalize(awnorm);
    }
    // It is enough to return awnorm and ret
    // sign(dist) = dot(grad, cell normal)
#if 0
    if (ret == 0)
    {
      for (int i = 0; i < 3; i++)
      {
        g[i] = awnorm[i];
      }
    }
#else
    // sign(dist) = dot(grad, cell normal)
    for (int i = 0; i < 3; i++)
    {
      g[i] = awnorm[i];
    }
#endif

    // Signed distance - we could skip this
    ret *= (vtkMath::Dot(g, awnorm) < 0.) ? 1. : -1.;

#if 0
    if (ret > 0.)
    {
      for (int i = 0; i < 3; i++)
      {
        g[i] = -g[i];
      }
    }
#endif
  }
  return intersectionFound;
}

int vtkImplicitPolyDataIntersection::IntersectWithLine(
  const double p1[3], const double p2[3], double& ret, double closestPoint[3], double g[3])
{
  vtkDataObject::AttributeTypes typ = vtkDataObject::NUMBER_OF_ATTRIBUTE_TYPES;
  return this->IntersectWithLine(p1, p2, ret, closestPoint, g, &typ);
}

//------------------------------------------------------------------------------
void vtkImplicitPolyDataIntersection::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os, indent);

  os << indent << "NoValue: " << this->NoValue << "\n";
  os << indent << "NoGradient: (" << this->NoGradient[0] << ", " << this->NoGradient[1] << ", "
     << this->NoGradient[2] << ")\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";

  if (this->Input)
  {
    os << indent << "Input : " << this->Input << "\n";
  }
  else
  {
    os << indent << "Input : (none)\n";
  }
}
