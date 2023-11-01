#pragma once

#include "config.h"

// For export macro
#if SPS_BUILD
#include "vtkSpsCommonDataModelModule.h"
#define VTKXXX_EXPORT VTKSPSCOMMONDATAMODEL_EXPORT
#else
#include "vtkICPModule.h"
#define VTKXXX_EXPORT VTKICP_EXPORT
#endif

#include "vtkLinearTransform.h"
#include "vtkPolyDataCorrespondenceFilter.h"

class vtkPolyData;
class vtkPolyDataCorrespondenceFilter;

class VTKXXX_EXPORT vtkICP : public vtkLinearTransform
{
public:
  static vtkICP* New();
  vtkTypeMacro(vtkICP, vtkLinearTransform)
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
    MetricPointToPlane = 1,
    MetricPlaneToPlane = 2, // TODO: Rename this to symmetrized
  };
  //@}

  //@{
  /**
   * Specify the source and target data sets.
   */
  void SetSource(vtkPolyData* source);
  void SetTarget(vtkPolyData* target);
  vtkGetObjectMacro(Source, vtkPolyData);
  vtkGetObjectMacro(Target, vtkPolyData);
  //@}

  //@{
  /**
   * Set/Get the maximum number of landmarks sampled in your dataset.
   * If your dataset is dense, then you will typically not need all the
   * points to compute the ICP transform. The default is 200.
   */
  vtkSetMacro(MaximumNumberOfLandmarks, int);
  vtkGetMacro(MaximumNumberOfLandmarks, int);
  //@}

  //@{
  /**
   * Set/Get the point-matching metric type
   */
  vtkSetClampMacro(Metric, int, MetricPointToPoint, MetricPlaneToPlane);
  vtkGetMacro(Metric, int);
  void SetMetricToPointToPoint() { this->SetMetric(MetricPointToPoint); }
  void SetMetricToPointToPlane() { this->SetMetric(MetricPointToPlane); }
  void SetMetricToPlaneToPlane() { this->SetMetric(MetricPlaneToPlane); }
  const char* GetMetricModeAsString();
  //@}

  //@{
  /**
   * Get the number of iterations since the last update
   */
  vtkGetMacro(NumberOfIterations, int);
  //@}

  //@{
  /**
   * Force the algorithm to check the mean distance between two iterations.
   * Default is Off.
   */
  vtkSetMacro(CheckMeanDistance, vtkTypeBool);
  vtkGetMacro(CheckMeanDistance, vtkTypeBool);
  vtkBooleanMacro(CheckMeanDistance, vtkTypeBool);
  //@}

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

  //  class VTKSPSCOMMONDATAMODEL_EXPORT vtkInternals;

protected:
  vtkICP();
  ~vtkICP() override;

  void InternalUpdate() override;

  /**
   * This method does no type checking, use DeepCopy instead.
   */
  void InternalDeepCopy(vtkAbstractTransform* transform) override;

  int MaximumNumberOfIterations;
  int NumberOfIterations;
  int MaximumNumberOfLandmarks;
  int Metric;
  vtkTypeBool CheckMeanDistance;
  vtkTypeBool StartByMatchingCentroids;
  double MaximumMeanDistance;
  double MeanDistance;
  double MaximumDistance;
  vtkPolyData* Source;
  vtkPolyData* Target;
  vtkPolyDataCorrespondenceFilter* Correspondences;

private:
  vtkICP(const vtkICP&) = delete;
  void operator=(const vtkICP&) = delete;
};
