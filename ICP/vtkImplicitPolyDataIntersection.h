/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImplicitPolyDataDistance.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkImplicitPolyDataDistance
 * @brief   Implicit function that computes the distance from a point x to the nearest point p on an
 * input vtkPolyData.
 *
 * Implicit function that computes the distance from a point x to the
 * nearest point p on an input vtkPolyData. The sign of the function
 * is set to the sign of the dot product between the angle-weighted
 * pseudonormal at the nearest surface point and the vector x - p.
 * Points interior to the geometry have a negative distance, points on
 * the exterior have a positive distance, and points on the input
 * vtkPolyData have a distance of zero. The gradient of the function
 * is the angle-weighted pseudonormal at the nearest point.
 *
 * Baerentzen, J. A. and Aanaes, H. (2005). Signed distance
 * computation using the angle weighted pseudonormal. IEEE
 * Transactions on Visualization and Computer Graphics, 11:243-253.
 *
 * This code was contributed in the VTK Journal paper:
 * "Boolean Operations on Surfaces in VTK Without External Libraries"
 * by Cory Quammen, Chris Weigle C., Russ Taylor
 * http://hdl.handle.net/10380/3262
 * http://www.midasjournal.org/browse/publication/797
 */

// TODO: Use different interface

#ifndef vtkImplicitPolyDataIntersection_h
#define vtkImplicitPolyDataIntersection_h

#include "vtkDataObject.h"
#include "vtkGenericCell.h" // For thread local storage
#include "vtkICPModule.h"   // For export macro
#include "vtkImplicitFunction.h"
#include "vtkObject.h"
#include "vtkSMPThreadLocalObject.h" // For thread local storage

class vtkCellLocator;
class vtkPolyData;

class VTKICP_EXPORT vtkImplicitPolyDataIntersection : public vtkObject
{
public:
  static vtkImplicitPolyDataIntersection* New();
  vtkTypeMacro(vtkImplicitPolyDataIntersection, vtkObject)

  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Return the MTime also considering the Input dependency.
   */
  vtkMTimeType GetMTime() override;

  /**
   * Set the input vtkPolyData used for the implicit function
   * evaluation.  Passes input through an internal instance of
   * vtkTriangleFilter to remove vertices and lines, leaving only
   * triangular polygons for evaluation as implicit planes.
   */
  void SetInput(vtkPolyData* input);

  ///@{
  /**
   * Set/get the function value to use if no input vtkPolyData
   * specified.
   */
  vtkSetMacro(NoValue, double);
  vtkGetMacro(NoValue, double);
  ///@}

  ///@{
  /**
   * Set/get the function gradient to use if no input vtkPolyData
   * specified.
   */
  vtkSetVector3Macro(NoGradient, double);
  vtkGetVector3Macro(NoGradient, double);
  ///@}

  ///@{
  /**
   * Set/get the closest point to use if no input vtkPolyData
   * specified.
   */
  vtkSetVector3Macro(NoClosestPoint, double);
  vtkGetVector3Macro(NoClosestPoint, double);
  ///@}

  ///@{
  /**
   * Set/get the tolerance used for the locator.
   */
  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);
  ///@}

  ///@{
  /**
   * Intersect the mesh with line from p1 to p2. Return the point and gradient for the
   * intersection as well as the distance t.
   */
  ///@}
  int IntersectWithLine(
    const double p1[3], const double p2[3], double& t, double x[3], double gradient[3]);

  int IntersectWithLine(const double p1[3], const double p2[3], double& t, double x[3],
    double gradient[3], vtkDataObject::AttributeTypes* type);

  int IntersectWithLine(
    const double p1[3], const double p2[3], double& t, double x[3], double gradient[3], int& id)
  {
    vtkDataObject::AttributeTypes type = vtkDataObject::NUMBER_OF_ATTRIBUTE_TYPES;
    int intersect = this->IntersectWithLine(p1, p2, t, x, gradient, &type);
    id = static_cast<int>(type);
    return intersect;
  }

protected:
  vtkImplicitPolyDataIntersection();
  ~vtkImplicitPolyDataIntersection() override;

  /**
   * Create default locator. Used to create one when none is specified.
   */
  void CreateDefaultLocator();

  double NoGradient[3];
  double NoClosestPoint[3];
  double NoValue;
  double Tolerance;

  vtkPolyData* Input;
  vtkCellLocator* Locator;
  vtkSMPThreadLocalObject<vtkGenericCell> TLCell;
  vtkSMPThreadLocalObject<vtkIdList> TLCellIds;

private:
  vtkImplicitPolyDataIntersection(const vtkImplicitPolyDataIntersection&) = delete;
  void operator=(const vtkImplicitPolyDataIntersection&) = delete;
};

#endif
