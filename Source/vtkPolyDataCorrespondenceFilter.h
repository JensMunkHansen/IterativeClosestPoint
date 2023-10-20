#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkICPModule.h" // For export macro

class vtkTransform;
class vtkImplicitPolyDataDistance2;
class vtkPolyDataNormals;

// TODO:
//  * Use vtkSmartPointer for the functor
//  * Support double precision

class VTKICP_EXPORT vtkPolyDataCorrespondenceFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPolyDataCorrespondenceFilter, vtkPolyDataAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkPolyDataCorrespondenceFilter* New();

  //@{
  /**
   *
   */
  enum SearchDirection
  {
    SourceToTarget = 0,
    TargetToSource = 1,
  };
  //@}

  vtkGetObjectMacro(Transform, vtkTransform);
  virtual void SetTransform(vtkTransform* transform);
  //@{
  /**
   * Set/Get the maximum number of landmarks sampled in your dataset.
   * If your dataset is dense, then you will typically not need all the
   * points to compute the ICP transform. The default is 200.
   */
  vtkSetMacro(MaximumNumberOfLandmarks, vtkIdType);
  vtkGetMacro(MaximumNumberOfLandmarks, vtkIdType);
  //@}

  vtkSetMacro(OutputPrecision, int);
  vtkGetMacro(OutputPrecision, int);

  //@{
  /**
   * Set/Get the maximum distance used for searching for correspondences.
   */
  vtkSetMacro(MaximumDistance, double);
  vtkGetMacro(MaximumDistance, double);
  //@}

  //@{
  /**
   * Search direction
   */
  vtkSetClampMacro(SearchDirection, int, SourceToTarget, TargetToSource);
  vtkGetMacro(SearchDirection, int);
  void SetSourceToTarget() { this->SetSearchDirection(SourceToTarget); }
  void SetTargetToSource() { this->SetSearchDirection(TargetToSource); }

  //@{
  /**
   *
   */
  vtkSetMacro(IncludeSourceNormals, vtkTypeBool);
  vtkGetMacro(IncludeSourceNormals, vtkTypeBool);
  //@}

protected:
  vtkPolyDataCorrespondenceFilter();
  ~vtkPolyDataCorrespondenceFilter() override;

  vtkMTimeType GetMTime() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  double MaximumDistance;
  vtkIdType MaximumNumberOfLandmarks;
  int SearchDirection;
  int OutputPrecision;
  vtkTypeBool IncludeSourceNormals;
  vtkTransform* Transform;
  vtkImplicitPolyDataDistance2* SourceLocator;
  vtkImplicitPolyDataDistance2* TargetLocator;
  vtkPolyDataNormals* SourceNormals;

private:
  vtkPolyDataCorrespondenceFilter(const vtkPolyDataCorrespondenceFilter&) = delete;
  void operator=(const vtkPolyDataCorrespondenceFilter) = delete;
};
