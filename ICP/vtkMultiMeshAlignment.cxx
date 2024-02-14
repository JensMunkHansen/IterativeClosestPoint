#include "vtkMultiMeshAlignment.h"
#include "icp/kabsch.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataCollection.h"
#include "vtkPolyDataCorrespondenceFilter.h"
#include "vtkTransform.h"
#include "vtkTransformCollection.h"

#include "multiMeshInternals.h"

#include "vtk_eigen.h"
#include VTK_EIGEN(Dense)

vtkStandardNewMacro(vtkMultiMeshAlignment);
vtkCxxSetObjectMacro(vtkMultiMeshAlignment, PolyDataCollection, vtkPolyDataCollection);
vtkCxxSetObjectMacro(vtkMultiMeshAlignment, TransformCollection, vtkTransformCollection);

//----------------------------------------------------------------------------
vtkMultiMeshAlignment::vtkMultiMeshAlignment()
  : vtkLinearTransform()
  , Internals(new vtkInternals(this))
{
  this->NumberOfIterations = 0;
  this->MaximumNumberOfIterations = 50;
  this->MeanDistance = 0.0;
  this->MaximumMeanDistance = 0.01;
  this->MaximumDistance = 10.0;
  this->Metric = vtkMultiMeshAlignment::MetricPointToPoint;
  this->Correspondences = vtkPolyDataCorrespondenceFilter::New();
  this->Correspondences->SetOutputPrecision(vtkAlgorithm::SINGLE_PRECISION);
  this->PolyDataCollection = vtkPolyDataCollection::New();
  this->TransformCollection = vtkTransformCollection::New();
}

//----------------------------------------------------------------------------
vtkMultiMeshAlignment::~vtkMultiMeshAlignment()
{
  this->PolyDataCollection->Delete();
  this->TransformCollection->Delete();
  this->Correspondences->Delete();
}

//----------------------------------------------------------------------------
vtkIdType vtkMultiMeshAlignment::GetNPolyData()
{
  return this->PolyDataCollection->GetNumberOfItems();
}

//----------------------------------------------------------------------------
bool vtkMultiMeshAlignment::PairPolyData(
  vtkPolyData* source, vtkPolyData* target, vtkTransform* initialTransform)
{
  int sourceIndex = this->PolyDataCollection->IndexOfFirstOccurence(source);
  int targetIndex = this->PolyDataCollection->IndexOfFirstOccurence(target);
  bool success = sourceIndex != -1 && targetIndex != -1;

  // 1. Use correspondence filter
  // 2. Criteria

  /*
        public static readonly PointSamplingParameters SamplingCorrespondingPoints = new
     PointSamplingParameters( step: 4, maxDist: 0.2f, minMatch: 5, interpolateNeighbours: true,
            minNormalDot: 0.95,
            useAbsoluteMetric: false);
  */

  // bite1 -> lower
  // bite1 -> upper
  // bite2 -> lower
  // bite2 -> upper

  return success;
}

//----------------------------------------------------------------------------
void vtkMultiMeshAlignment::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);
  if (this->PolyDataCollection)
  {
    this->PolyDataCollection->InitTraversal();
    vtkPolyData* polyData = this->PolyDataCollection->GetNextItem();
    vtkIdType polyIndex = 0;
    vtkIdType transformIndex = 0;
    this->TransformCollection->InitTraversal();
    vtkTransform* transform = this->TransformCollection->GetNextItem();
    while (polyData != nullptr)
    {
      os << indent << "PolyData[" << polyIndex << "]:\n";
      os << indent << indent << *polyData;
      polyData = this->PolyDataCollection->GetNextItem();
      if (transform != nullptr)
      {
        os << indent << "Transform[" << transformIndex << "]:\n";
        os << indent << indent << *transform;
      }
      polyIndex++;
      transformIndex++;
    }
  }
  else
  {
    os << indent << "PolyData: (none)\n";
  }

  os << indent << "MaximumNumberOfIterations: " << this->MaximumNumberOfIterations << "\n";
  os << indent << "MaximumMeanDistance: " << this->MaximumMeanDistance << "\n";
  os << indent << "MaximumDistance: " << this->MaximumDistance << "\n";
  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
  os << indent << "MeanDistance: " << this->MeanDistance << "\n";
  os << indent << "Metric: " << this->GetMetricModeAsString() << "\n";
}

//----------------------------------------------------------------------------
const char* vtkMultiMeshAlignment::GetMetricModeAsString()
{
  const char* retval = nullptr;

  switch (this->Metric)
  {
    case MetricPointToPoint:
    default:
      retval = "Point-to-point";
      break;
  }
  return retval;
}

//----------------------------------------------------------------------------
void vtkMultiMeshAlignment::Inverse()
{
  /*
  vtkPolyData* tmp1 = this->Source;
  this->Source = this->Target;
  this->Target = tmp1;
  this->Modified();
  */
}

//----------------------------------------------------------------------------
void vtkMultiMeshAlignment::InternalDeepCopy(vtkAbstractTransform* transform)
{
  Superclass::InternalDeepCopy(transform);

  vtkMultiMeshAlignment* t = (vtkMultiMeshAlignment*)transform;

  this->SetMaximumNumberOfIterations(t->GetMaximumNumberOfIterations());
  this->SetMaximumMeanDistance(t->GetMaximumMeanDistance());
  this->SetMaximumNumberOfLandmarks(t->GetMaximumNumberOfLandmarks());
  this->SetMaximumDistance(t->GetMaximumDistance());

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkMultiMeshAlignment::InternalUpdate()
{
  vtkDebugMacro(<< __FUNCTION__);
  if (this->PolyDataCollection->GetNumberOfItems() == this->TransformCollection->GetNumberOfItems())
  {
    // We can update - expensive
  }
}

bool vtkMultiMeshAlignment::TestKabsch()
{
  icp::Kabsch<float>* kabsch = new icp::Kabsch<float>();
  size_t nPoints = 10;
  std::vector<vtkeigen::Vector3f> target(nPoints);
  std::vector<vtkeigen::Vector3f> source(nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    target[i] = vtkeigen::Vector3f(1, 2, 3);
    source[i] = vtkeigen::Vector3f(1, 2, 3);
  }
  vtkeigen::Affine3f tx;
  float residuals = 0;
  kabsch->Update(target, source, tx, residuals);
  delete kabsch;
  cout << tx.affine();
  return true;
}

//----------------------------------------------------------------------------
vtkAbstractTransform* vtkMultiMeshAlignment::MakeTransform()
{
  return vtkMultiMeshAlignment::New();
}
