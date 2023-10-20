#pragma once

#include "vtkPolyDataCorrespondenceFilter.h"
#include "vtkICPModule.h" // For export macro

#include "vtkLinearTransform.h"

class vtkPolyData;
class vtkPolyDataCorrespondenceFilter;

// TODO: Handle general vtkPointSet (not just vtkPolyData)

class VTKICP_EXPORT vtkICP : public vtkLinearTransform
{
public:
  static vtkICP* New();
  vtkTypeMacro(vtkICP, vtkLinearTransform)
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(MaximumNumberOfIterations, int);
  vtkSetMacro(MaximumNumberOfIterations, int);

  enum MetricType : int
  {
    MetricPointToPoint = 0,
    MetricPointToPlane = 1,
  };
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
   * Specify the mean distance mode. This mode expresses how the mean
   * distance is computed. The RMS mode is the square root of the average
   * of the sum of squares of the closest point distances. The Absolute
   * Value mode is the mean of the sum of absolute values of the closest
   * point distances. The default is VTK_ICP_MODE_RMS
   */
  vtkSetClampMacro(Metric, int, MetricPointToPoint, MetricPointToPlane);
  vtkGetMacro(Metric, int);
  void SetMetricToPointToPoint() { this->SetMetric(MetricPointToPoint); }
  void SetMetricToPointToPlane() { this->SetMetric(MetricPointToPlane); }
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
   * Starts the process by translating source centroid to target centroid.
   * The default is Off.
   */
  vtkSetMacro(StartByMatchingCentroids, vtkTypeBool);
  vtkGetMacro(StartByMatchingCentroids, vtkTypeBool);
  vtkBooleanMacro(StartByMatchingCentroids, vtkTypeBool);
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
  vtkPolyData* Source;
  vtkPolyData* Target;
  vtkPolyDataCorrespondenceFilter* Correspondences;

private:
  vtkICP(const vtkICP&) = delete;
  void operator=(const vtkICP&) = delete;
};
