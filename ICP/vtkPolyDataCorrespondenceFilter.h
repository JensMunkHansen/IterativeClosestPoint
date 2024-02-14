#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkICPModule.h" // For export macro

class vtkTransform;
class vtkImplicitPolyDataDistance2;
class vtkPolyDataNormals;

class VTKICP_EXPORT vtkPolyDataCorrespondenceFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPolyDataCorrespondenceFilter, vtkPolyDataAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkPolyDataCorrespondenceFilter* New();

  //@{
  /**
   * Search direction for correspondences. Source is input port 0,
   * target is input port 1.
   */
  enum SearchDirection
  {
    SourceToTarget = 0,
    TargetToSource = 1,
  };
  //@}

  //@{
  /**
   * Set/Get an implicit transform applied to the source from input
   * port 0. This transform is applied before correspondences are
   * found.
   */
  vtkGetObjectMacro(Transform, vtkTransform);
  virtual void SetTransform(vtkTransform* transform);
  //@}

  //@{
  /**
   * Set/Get the maximum number of landmarks sampled in your dataset.
   * If your dataset is dense, then you will typically not need all the
   * points to compute the ICP transform. The default is VTK_INT_MAX
   */
  vtkSetMacro(MaximumNumberOfLandmarks, vtkIdType);
  vtkGetMacro(MaximumNumberOfLandmarks, vtkIdType);
  //@}

  //@{
  /**
   * Output precision for the points and normals. The default is SINGLE_PRECISION
   */
  vtkSetMacro(OutputPrecision, int);
  vtkGetMacro(OutputPrecision, int);
  //@}

  //@{
  /**
   * Set/Get the maximum distance used for searching for correspondences. Default is 10.0
   */
  vtkSetMacro(MaximumDistance, double);
  vtkGetMacro(MaximumDistance, double);
  //@}

  //@{
  /**
   * Set/Get the minimum value for the dot product of normals. Default is 0.0
   */
  vtkSetMacro(MinNormalDot, double);
  vtkGetMacro(MinNormalDot, double);
  //@}

  //@{
  /**
   * Search direction. The default is from source (port 0) to target (port 1)
   */
  vtkSetClampMacro(SearchDirection, int, SourceToTarget, TargetToSource);
  vtkGetMacro(SearchDirection, int);
  void SetSourceToTarget() { this->SetSearchDirection(SourceToTarget); }
  void SetTargetToSource() { this->SetSearchDirection(TargetToSource); }

  //@{
  /**
   * Include origin normals.
   */
  vtkSetMacro(IncludeOriginNormals, vtkTypeBool);
  vtkGetMacro(IncludeOriginNormals, vtkTypeBool);
  //@}

protected:
  vtkPolyDataCorrespondenceFilter();
  ~vtkPolyDataCorrespondenceFilter() override;

  vtkMTimeType GetMTime() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  double MaximumDistance;
  double MinNormalDot;
  vtkIdType MaximumNumberOfLandmarks;
  int SearchDirection;
  int OutputPrecision;
  vtkTypeBool IncludeOriginNormals;
  vtkTransform* Transform;
  vtkImplicitPolyDataDistance2* SourceLocator;
  vtkImplicitPolyDataDistance2* TargetLocator;
  vtkPolyDataNormals* SourceNormals;

private:
  vtkPolyDataCorrespondenceFilter(const vtkPolyDataCorrespondenceFilter&) = delete;
  void operator=(const vtkPolyDataCorrespondenceFilter) = delete;
};
