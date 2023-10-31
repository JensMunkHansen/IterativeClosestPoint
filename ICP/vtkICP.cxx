#include "vtkICP.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataCorrespondenceFilter.h"
#include "vtkTransform.h"

#include "vtk_eigen.h"
#include VTK_EIGEN(Dense)

#include <icp/point2plane.hpp>
#include <icp/point2point.hpp>
#include <icp/plane2plane.hpp>

using vtkeigen::Map;

vtkStandardNewMacro(vtkICP);

vtkCxxSetObjectMacro(vtkICP, Source, vtkPolyData);
vtkCxxSetObjectMacro(vtkICP, Target, vtkPolyData);

//----------------------------------------------------------------------------
vtkICP::vtkICP()
  : vtkLinearTransform()
{
  this->NumberOfIterations = 0;
  this->MaximumNumberOfIterations = 50;
  this->CheckMeanDistance = false;
  this->MeanDistance = 0.0;
  this->MaximumMeanDistance = 0.01;
  this->MaximumDistance = 10.0;
  this->Metric = vtkICP::MetricPointToPlane;
  this->Correspondences = vtkPolyDataCorrespondenceFilter::New();
  this->Correspondences->SetOutputPrecision(vtkAlgorithm::SINGLE_PRECISION);
  this->Source = this->Target = nullptr;
}

//----------------------------------------------------------------------------
vtkICP::~vtkICP()
{
  this->Correspondences->Delete();
  if (this->Source)
    Source->UnRegister(this);
  if (this->Target)
    Target->UnRegister(this);
}

//----------------------------------------------------------------------------
void vtkICP::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);
  if (this->Source)
  {
    os << indent << "Source: " << this->Source << "\n";
  }
  else
  {
    os << indent << "Source: (none)\n";
  }

  if (this->Target)
  {
    os << indent << "Target: " << this->Target << "\n";
  }
  else
  {
    os << indent << "Target: (none)\n";
  }

  os << indent << "MaximumNumberOfIterations: " << this->MaximumNumberOfIterations << "\n";
  os << indent << "CheckMeanDistance: " << this->CheckMeanDistance << "\n";
  os << indent << "MaximumMeanDistance: " << this->MaximumMeanDistance << "\n";
  os << indent << "MaximumDistance: " << this->MaximumDistance << "\n";
  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
  os << indent << "MeanDistance: " << this->MeanDistance << "\n";
  os << indent << "Metric: " << this->GetMetricModeAsString() << "\n";
}

//----------------------------------------------------------------------------
const char* vtkICP::GetMetricModeAsString()
{
  if (this->Metric == MetricPointToPoint)
  {
    return "Point-to-point";
  }
  else if (this->Metric == MetricPointToPlane)
  {
    return "Point-to-plane";
  }
  else
  {
    return "Plane-to-plane";
  }
}

//----------------------------------------------------------------------------
void vtkICP::Inverse()
{
  vtkPolyData* tmp1 = this->Source;
  this->Source = this->Target;
  this->Target = tmp1;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkICP::InternalDeepCopy(vtkAbstractTransform* transform)
{
  Superclass::InternalDeepCopy(transform);

  vtkICP* t = (vtkICP*)transform;

  this->SetSource(t->GetSource());
  this->SetTarget(t->GetTarget());
  this->SetMaximumNumberOfIterations(t->GetMaximumNumberOfIterations());
  this->SetCheckMeanDistance(t->GetCheckMeanDistance());
  this->SetMaximumMeanDistance(t->GetMaximumMeanDistance());
  this->SetMaximumNumberOfLandmarks(t->GetMaximumNumberOfLandmarks());

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkICP::InternalUpdate()
{
  vtkDebugMacro(<< __FUNCTION__);
  if (this->Source == nullptr || !this->Source->GetNumberOfPoints())
  {
    vtkErrorMacro(<< "Can't execute with nullptr or empty input");
    return;
  }

  if (this->Target == nullptr || !this->Target->GetNumberOfPoints())
  {
    vtkErrorMacro(<< "Can't execute with nullptr or empty target");
    return;
  }

  if (this->Target->GetPoints()->GetDataType() != VTK_FLOAT ||
    this->Source->GetPoints()->GetDataType() != VTK_FLOAT)
  {
    vtkErrorMacro(<< "We only support single-precision points for now");
    return;
  }

  this->Correspondences->SetInputData(0, this->Source);
  this->Correspondences->SetInputData(1, this->Target);
  this->Correspondences->SetMaximumNumberOfLandmarks(this->GetMaximumNumberOfLandmarks());
  this->Correspondences->SetMaximumDistance(this->GetMaximumDistance());
  
  icp::ICP<float, float, 3>* icp = nullptr;
  if (this->Metric == MetricPointToPlane)
    icp = new icp::PointToPlane<float, float, 3>;
  else if (this->Metric == MetricPointToPoint)
    icp = new icp::PointToPoint<float, float, 3>;
  else
  {
    icp = new icp::PlaneToPlane<float, float, 3>;
    this->Correspondences->SetIncludeOriginNormals(true);
  }

  vtkNew<vtkTransform> transform;
  transform->PostMultiply();
  this->Correspondences->SetTransform(transform);

  this->NumberOfIterations = 0;
  vtkFloatArray* sourcePointDataArray = nullptr;
  vtkFloatArray* sourceNormalsDataArray = nullptr;
  vtkFloatArray* targetPointDataArray = nullptr;
  vtkFloatArray* targetNormalsDataArray = nullptr;
  vtkDoubleArray* targetNormalsDataArrayDouble = nullptr;

  if (this->Metric == MetricPlaneToPlane)
  {
    // TODO: Add check for source normals
  }

  
  for (vtkIdType it = 0; it < this->MaximumNumberOfIterations; it++)
  {
    this->Correspondences->Update();

    sourcePointDataArray =
      vtkFloatArray::SafeDownCast(this->Correspondences->GetOutput(0)->GetPoints()->GetData());
    sourceNormalsDataArray = vtkFloatArray::SafeDownCast(
        this->Correspondences->GetOutput(0)->GetPointData()->GetArray(
            "Normals"));
    targetPointDataArray =
      vtkFloatArray::SafeDownCast(this->Correspondences->GetOutput(1)->GetPoints()->GetData());
    targetNormalsDataArray = vtkFloatArray::SafeDownCast(
      this->Correspondences->GetOutput(1)->GetPointData()->GetArray("Normals"));
    targetNormalsDataArrayDouble = vtkDoubleArray::SafeDownCast(
      this->Correspondences->GetOutput(1)->GetPointData()->GetArray("Normals"));

    vtkeigen::Matrix3Xf source = Map<vtkeigen::Matrix3Xf>(
      sourcePointDataArray->GetPointer(0), 3, sourcePointDataArray->GetNumberOfTuples());
    vtkeigen::Matrix3Xf target = Map<vtkeigen::Matrix3Xf>(
      targetPointDataArray->GetPointer(0), 3, targetPointDataArray->GetNumberOfTuples());
    vtkeigen::Matrix3Xf normals = Map<vtkeigen::Matrix3Xf>(
      targetNormalsDataArray->GetPointer(0), 3, targetNormalsDataArray->GetNumberOfTuples());

    vtkeigen::Matrix3Xf sourceNormals;
    
    if (this->Metric == MetricPlaneToPlane)
      sourceNormals = Map<vtkeigen::Matrix3Xf>(
        sourceNormalsDataArray->GetPointer(0), 3, sourceNormalsDataArray->GetNumberOfTuples());
    
    vtkeigen::VectorXf W = vtkeigen::VectorXf::Ones(source.cols());
    vtkeigen::VectorXf U = vtkeigen::VectorXf::Zero(source.cols());

    // Here we can check mean distance
    if (this->CheckMeanDistance)
    {
      this->MeanDistance = (source - target).colwise().norm().mean();
      if (this->MeanDistance <= this->MaximumMeanDistance)
      {
        break;
      }
    }
    auto affine = icp->Update(source, target, sourceNormals, normals, W);
    vtkNew<vtkMatrix4x4> mat;
    auto eigenMat = affine.matrix();
    this->NumberOfIterations++;

    // We need to transpose anyway
    for (size_t c = 0; c < eigenMat.cols(); c++)
      for (size_t r = 0; r < eigenMat.rows(); r++)
        mat->SetElement(r, c, eigenMat(r, c));

    this->Correspondences->GetTransform()->Concatenate(mat);
    this->Correspondences->GetTransform()->Update();
  }

  delete icp;

  this->Matrix->DeepCopy(this->Correspondences->GetTransform()->GetMatrix());
}

//----------------------------------------------------------------------------
vtkAbstractTransform* vtkICP::MakeTransform()
{
  return vtkICP::New();
}
