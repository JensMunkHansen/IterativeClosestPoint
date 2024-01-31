#pragma once

#include "vtkICPModule.h"
#include "vtkLinearTransform.h"

class vtkPolyData;
class vtkPolyDataCollection;
class vtkTransformCollection;
class vtkPolyDataCorrespondenceFilter;
class vtkTransform;

class VTKICP_EXPORT vtkMultiMeshAlignment : public vtkLinearTransform
{
public:
  static vtkMultiMeshAlignment* New();
  vtkTypeMacro(vtkMultiMeshAlignment, vtkLinearTransform)
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(MaximumNumberOfIterations, int);
  vtkSetMacro(MaximumNumberOfIterations, int);

  //@{
  /**
   * Point-matching metric types
   */
  enum MetricType : int
  {
    MetricPointToPoint = 0,
  };
  //@}

  //@{
  /**
   * Specify the source and target data sets.
   */
  void SetPolyDataCollection(vtkPolyDataCollection* collection);
  //@}

  //@{
  /**
   * Specify the source and target data sets.
   */
  void SetTransformCollection(vtkTransformCollection* collection);
  //@}

  //@{
  /**
   * Pair polydata, we register source to target
   */
  bool PairPolyData(vtkPolyData* source, vtkPolyData* target, vtkTransform* initialTransform);
  //@}

  //@{
  /**
   * Set/Get the maximum number of landmarks sampled in your dataset.
   * If your dataset is dense, then you will typically not need all the
   * points to compute the alignment. The default is 200.
   */
  vtkSetMacro(MaximumNumberOfLandmarks, int);
  vtkGetMacro(MaximumNumberOfLandmarks, int);
  //@}

  //@{
  /**
   * Set/Get the point-matching metric type
   */
  vtkSetClampMacro(Metric, int, MetricPointToPoint, MetricPointToPoint);
  vtkGetMacro(Metric, int);
  void SetMetricToPointToPoint() { this->SetMetric(MetricPointToPoint); }
  const char* GetMetricModeAsString();
  //@}

  //@{
  /**
   * Get the number of iterations since the last update
   */
  vtkGetMacro(NumberOfIterations, int);
  //@}

  vtkIdType GetNPolyData();

  /**
   * Invert the transformation.  This is done by switching the
   * source and target.
   */
  void Inverse() override;

  //@{
  /**
   * Set/Get the maximum mean distance between two iteration. If the mean
   * distance is lower than this, the convergence stops. The default
   * is 0.01.
   */
  vtkSetMacro(MaximumMeanDistance, double);
  vtkGetMacro(MaximumMeanDistance, double);
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
   * Get the mean distance between the last two iterations.
   */
  vtkGetMacro(MeanDistance, double);
  //@}

  /**
   * Make another transform of the same type.
   */
  vtkAbstractTransform* MakeTransform() override;

  bool TestKabsch();

  // PIMPL
  class VTKICP_EXPORT vtkInternals;

protected:
  vtkMultiMeshAlignment();
  ~vtkMultiMeshAlignment() override;

  void InternalUpdate() override;

  /**
   * This method does no type checking, use DeepCopy instead.
   */
  void InternalDeepCopy(vtkAbstractTransform* transform) override;

  int MaximumNumberOfIterations;
  int NumberOfIterations;
  int MaximumNumberOfLandmarks;
  int Metric;
  double MaximumMeanDistance;
  double MeanDistance;
  double MaximumDistance;
  vtkPolyDataCollection* PolyDataCollection;
  vtkTransformCollection* TransformCollection;
  vtkPolyDataCorrespondenceFilter* Correspondences;

  vtkInternals* Internals;
  friend class vtkInternals;

private:
  vtkMultiMeshAlignment(const vtkMultiMeshAlignment&) = delete;
  void operator=(const vtkMultiMeshAlignment&) = delete;
};
