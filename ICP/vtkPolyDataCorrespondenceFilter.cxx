#include "vtkPolyDataCorrespondenceFilter.h"

#include "vtkImplicitPolyDataDistance2.h"
#include "vtkObjectFactory.h"
#include "vtkSMPThreadLocalObject.h"
#include "vtkSMPTools.h"
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTransform.h>

namespace
{
template <typename TPointsArrayIn, typename TPointsArrayOut, typename TVectorArrayOut>
struct Functor
{
  struct vtkLocalDataType
  {
    TPointsArrayOut* Origins;
    TPointsArrayOut* Landmarks;
    TVectorArrayOut* Normals;
    int LastIndex;
    vtkLocalDataType()
      : Origins(nullptr)
      , Landmarks(nullptr)
      , Normals(nullptr)
    {
    }
  };

  vtkSMPThreadLocalObject<TPointsArrayOut> TLOrigins;
  vtkSMPThreadLocalObject<TPointsArrayOut> TLLandmarks;
  vtkSMPThreadLocalObject<TVectorArrayOut> TLNormals;
  vtkSMPThreadLocal<int> TLLastIndex;
  vtkSMPThreadLocal<vtkLocalDataType> LocalData;

  TPointsArrayIn* Points;
  TPointsArrayOut* Origins;
  TPointsArrayOut* Landmarks;
  TVectorArrayOut* Normals;
  vtkIdType Step;
  // TODO: Consider exposing the active source and target locator on filter
  vtkImplicitPolyDataDistance2* Locator;
  vtkPolyDataCorrespondenceFilter* Filter;

  Functor(TPointsArrayIn* points, TPointsArrayOut* origins, TPointsArrayOut* landmarks,
    TVectorArrayOut* normals, vtkIdType step, vtkImplicitPolyDataDistance2* locator,
    vtkPolyDataCorrespondenceFilter* filter)
    : Points(points)
    , Origins(origins)
    , Landmarks(landmarks)
    , Normals(normals)
    , Step(step)
    , Locator(locator)
    , Filter(filter)
  {
  }
  void Initialize()
  {
    auto& localData = this->LocalData.Local();
    auto& newOrigins = this->TLOrigins.Local();
    auto& newLandmarks = this->TLLandmarks.Local();
    auto& newNormals = this->TLNormals.Local();
    int& newLastIndex = this->TLLastIndex.Local();
    newLastIndex = 0;

    // The maximum number of correspondences for all threads
    vtkIdType estSize = std::ceil((float)this->Points->GetNumberOfTuples() / this->Step);

    // Allocate per thread arrays
    newOrigins->SetNumberOfComponents(3);
    newLandmarks->SetNumberOfComponents(3);
    newNormals->SetNumberOfComponents(3);
    newOrigins->SetNumberOfTuples(estSize);
    newLandmarks->SetNumberOfTuples(estSize);
    newNormals->SetNumberOfTuples(estSize);

    localData.Landmarks = newLandmarks;
    localData.Normals = newNormals;
    localData.Origins = newOrigins;
    localData.LastIndex = newLastIndex;
  }

  void operator()(vtkIdType begin, vtkIdType end)
  {

    auto& localData = this->LocalData.Local();

    using VecTypeOut = typename vtkDataArrayAccessor<TVectorArrayOut>::APIType;
    using PointTypeOut = typename vtkDataArrayAccessor<TPointsArrayOut>::APIType;

    auto landmarks = vtk::DataArrayTupleRange<3>(localData.Landmarks);
    auto origins = vtk::DataArrayTupleRange<3>(localData.Origins);
    auto normals = vtk::DataArrayTupleRange<3>(localData.Normals);

    vtkDataArrayAccessor<TPointsArrayIn> points(this->Points);

    const int comps = localData.Normals->GetNumberOfComponents();
    int& lastIndex = localData.LastIndex;

    double signedDistance;
    double originPoint[3];
    double closestPoint[3];
    double normal[3];
    vtkDataObject::AttributeTypes type = vtkDataObject::NUMBER_OF_ATTRIBUTE_TYPES;

    for (; begin < end; ++begin)
    {
      if (begin % Step != 0)
        continue;

      // vtkImplicitPolyDataDistance uses double precision
      originPoint[0] = static_cast<double>(points.Get(begin, 0));
      originPoint[1] = static_cast<double>(points.Get(begin, 1));
      originPoint[2] = static_cast<double>(points.Get(begin, 2));
      Filter->GetTransform()->InternalTransformPoint(originPoint, originPoint);
      signedDistance = Locator->ClosestPointAndNormal(originPoint, normal, closestPoint);

      if (std::fabs(signedDistance) < Filter->GetMaximumDistance())
      {
        auto landmark = landmarks[lastIndex];
        auto origin = origins[lastIndex];
        auto normal = normals[lastIndex];

        for (int c = 0; c < comps; c++)
        {
          origin[c] = static_cast<PointTypeOut>(originPoint[c]);
          landmark[c] = static_cast<PointTypeOut>(closestPoint[c]);
          normal[c] = static_cast<VecTypeOut>(normal[c]);
        }
        lastIndex++;
      }
    }
  }
  void Reduce()
  {
    using TLSIter = typename vtkSMPThreadLocal<vtkLocalDataType>::iterator;
    TLSIter end = this->LocalData.end();
    vtkIdType dstStart = 0;

    vtkIdType nTotalCount = 0;
    for (TLSIter itr = this->LocalData.begin(); itr != end; ++itr)
    {
      nTotalCount += (*itr).LastIndex;
    }
    // Allocate outputs
    this->Landmarks->SetNumberOfComponents(3);
    this->Landmarks->SetNumberOfTuples(nTotalCount);
    this->Origins->SetNumberOfComponents(3);
    this->Origins->SetNumberOfTuples(nTotalCount);
    this->Normals->SetNumberOfComponents(3);
    this->Normals->SetNumberOfTuples(nTotalCount);

    for (TLSIter itr = this->LocalData.begin(); itr != end; ++itr)
    {
      auto origins = (*itr).Origins;
      auto landmarks = (*itr).Landmarks;
      auto normals = (*itr).Normals;
      vtkIdType nPoints = (*itr).LastIndex;
      this->Origins->InsertTuples(dstStart, nPoints, 0, origins);
      this->Landmarks->InsertTuples(dstStart, nPoints, 0, landmarks);
      this->Normals->InsertTuples(dstStart, nPoints, 0, normals);
      dstStart += nPoints;
    }
  }
};

template <typename TPointsArrayIn, typename TPointsArrayOut, typename TVectorArrayIn,
  typename TVectorArrayOut>
struct NormalsFunctor
{
  struct vtkLocalDataType
  {
    TPointsArrayOut* Origins;
    TVectorArrayOut* OriginNormals;
    TPointsArrayOut* Landmarks;
    TVectorArrayOut* LandmarkNormals;
    int LastIndex;
    vtkLocalDataType()
      : Origins(nullptr)
      , OriginNormals(nullptr)
      , Landmarks(nullptr)
      , LandmarkNormals(nullptr)
    {
    }
  };

  vtkSMPThreadLocalObject<TPointsArrayOut> TLOrigins;
  vtkSMPThreadLocalObject<TPointsArrayOut> TLOriginNormals;
  vtkSMPThreadLocalObject<TPointsArrayOut> TLLandmarks;
  vtkSMPThreadLocalObject<TVectorArrayOut> TLLandmarkNormals;
  vtkSMPThreadLocal<int> TLLastIndex;
  vtkSMPThreadLocal<vtkLocalDataType> LocalData;

  TPointsArrayIn* Points;
  TVectorArrayIn* PointNormals;
  TPointsArrayOut* Origins;
  TVectorArrayOut* OriginNormals;
  TPointsArrayOut* Landmarks;
  TVectorArrayOut* LandmarkNormals;
  vtkIdType Step;
  vtkImplicitPolyDataDistance2* Locator;
  vtkPolyDataCorrespondenceFilter* Filter;

  NormalsFunctor(TPointsArrayIn* points, TVectorArrayIn* pointNormals, TPointsArrayOut* origins,
    TVectorArrayOut* originNormals, TPointsArrayOut* landmarks, TVectorArrayOut* landmarkNormals,
    vtkIdType step, vtkImplicitPolyDataDistance2* locator, vtkPolyDataCorrespondenceFilter* filter)
    : Points(points)
    , PointNormals(pointNormals)
    , Origins(origins)
    , OriginNormals(originNormals)
    , Landmarks(landmarks)
    , LandmarkNormals(landmarkNormals)
    , Step(step)
    , Locator(locator)
    , Filter(filter)
  {
  }
  void Initialize()
  {
    auto& localData = this->LocalData.Local();
    auto& newOrigins = this->TLOrigins.Local();
    auto& newOriginNormals = this->TLOriginNormals.Local();
    auto& newLandmarks = this->TLLandmarks.Local();
    auto& newLandmarkNormals = this->TLLandmarkNormals.Local();
    int& newLastIndex = this->TLLastIndex.Local();
    newLastIndex = 0;

    // The maximum number of correspondences for all threads
    vtkIdType estSize = std::ceil((float)this->Points->GetNumberOfTuples() / this->Step);

    // Allocate per thread arrays
    newOrigins->SetNumberOfComponents(3);
    newOriginNormals->SetNumberOfComponents(3);
    newLandmarks->SetNumberOfComponents(3);
    newLandmarkNormals->SetNumberOfComponents(3);
    newOrigins->SetNumberOfTuples(estSize);
    newOriginNormals->SetNumberOfTuples(estSize);
    newLandmarks->SetNumberOfTuples(estSize);
    newLandmarkNormals->SetNumberOfTuples(estSize);

    localData.Landmarks = newLandmarks;
    localData.LandmarkNormals = newLandmarkNormals;
    localData.Origins = newOrigins;
    localData.OriginNormals = newOriginNormals;
    localData.LastIndex = newLastIndex;
  }

  void operator()(vtkIdType begin, vtkIdType end)
  {
    auto& localData = this->LocalData.Local();

    using VecTypeOut = typename vtkDataArrayAccessor<TVectorArrayOut>::APIType;
    using PointTypeOut = typename vtkDataArrayAccessor<TPointsArrayOut>::APIType;

    auto origins = vtk::DataArrayTupleRange<3>(localData.Origins);
    auto originNormals = vtk::DataArrayTupleRange<3>(localData.OriginNormals);
    auto landmarks = vtk::DataArrayTupleRange<3>(localData.Landmarks);
    auto landmarkNormals = vtk::DataArrayTupleRange<3>(localData.LandmarkNormals);
    vtkDataArrayAccessor<TPointsArrayIn> points(this->Points);
    vtkDataArrayAccessor<TPointsArrayIn> pointNormals(this->PointNormals);

    const int comps = localData.LandmarkNormals->GetNumberOfComponents();
    int& lastIndex = localData.LastIndex;

    double signedDistance;
    double originPoint[3];
    double closestPoint[3];
    double normal[3];

    for (; begin < end; ++begin)
    {
      if (begin % Step != 0)
        continue;

      // vtkImplicitPolyDataDistance uses double precision
      originPoint[0] = static_cast<double>(points.Get(begin, 0));
      originPoint[1] = static_cast<double>(points.Get(begin, 1));
      originPoint[2] = static_cast<double>(points.Get(begin, 2));

      Filter->GetTransform()->InternalTransformPoint(originPoint, originPoint);
      signedDistance = Locator->ClosestPointAndNormal(originPoint, normal, closestPoint);

      if (std::fabs(signedDistance) < Filter->GetMaximumDistance())
      {
        auto origin = origins[lastIndex];
        auto originNormal = originNormals[lastIndex];
        auto landmark = landmarks[lastIndex];
        auto landmarkNormal = landmarkNormals[lastIndex];

        for (int c = 0; c < comps; c++)
        {
          origin[c] = static_cast<PointTypeOut>(originPoint[c]);
          originNormal[c] = static_cast<VecTypeOut>(pointNormals.Get(begin, c));
          landmark[c] = static_cast<PointTypeOut>(closestPoint[c]);
          landmarkNormal[c] = static_cast<VecTypeOut>(normal[c]);
        }
        lastIndex++;
      }
    }
  }
  void Reduce()
  {
    using TLSIter = typename vtkSMPThreadLocal<vtkLocalDataType>::iterator;
    TLSIter end = this->LocalData.end();
    vtkIdType dstStart = 0;

    vtkIdType nTotalCount = 0;
    for (TLSIter itr = this->LocalData.begin(); itr != end; ++itr)
    {
      nTotalCount += (*itr).LastIndex;
    }
    // Allocate outputs
    this->Landmarks->SetNumberOfComponents(3);
    this->Landmarks->SetNumberOfTuples(nTotalCount);
    this->Origins->SetNumberOfComponents(3);
    this->Origins->SetNumberOfTuples(nTotalCount);
    this->OriginNormals->SetNumberOfComponents(3);
    this->OriginNormals->SetNumberOfTuples(nTotalCount);
    this->LandmarkNormals->SetNumberOfComponents(3);
    this->LandmarkNormals->SetNumberOfTuples(nTotalCount);

    for (TLSIter itr = this->LocalData.begin(); itr != end; ++itr)
    {
      auto origins = (*itr).Origins;
      auto originNormals = (*itr).OriginNormals;
      auto landmarks = (*itr).Landmarks;
      auto landmarkNormals = (*itr).LandmarkNormals;
      vtkIdType nPoints = (*itr).LastIndex;
      this->Origins->InsertTuples(dstStart, nPoints, 0, origins);
      this->OriginNormals->InsertTuples(dstStart, nPoints, 0, originNormals);
      this->Landmarks->InsertTuples(dstStart, nPoints, 0, landmarks);
      this->LandmarkNormals->InsertTuples(dstStart, nPoints, 0, landmarkNormals);
      dstStart += nPoints;
    }
  }
};
}

vtkStandardNewMacro(vtkPolyDataCorrespondenceFilter);

vtkCxxSetObjectMacro(vtkPolyDataCorrespondenceFilter, Transform, vtkTransform);

vtkPolyDataCorrespondenceFilter::vtkPolyDataCorrespondenceFilter()
{
  this->MaximumDistance = 10.0;
  this->MaximumNumberOfLandmarks = VTK_INT_MAX;
  this->OutputPrecision = vtkAlgorithm::DEFAULT_PRECISION;
  this->SourceLocator = vtkImplicitPolyDataDistance2::New();
  this->TargetLocator = vtkImplicitPolyDataDistance2::New();
  this->SourceNormals = vtkPolyDataNormals::New();
  this->SourceNormals->SplittingOff();
  this->SourceNormals->ConsistencyOn();
  this->SourceNormals->ComputePointNormalsOn();

  this->IncludeOriginNormals = false;
  this->SearchDirection = SourceToTarget;
  this->Transform = vtkTransform::New();
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

vtkPolyDataCorrespondenceFilter::~vtkPolyDataCorrespondenceFilter()
{
  this->SourceLocator->Delete();
  this->TargetLocator->Delete();
  this->SourceNormals->Delete();
  if (this->Transform)
    this->Transform->Delete();
}

void vtkPolyDataCorrespondenceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

vtkMTimeType vtkPolyDataCorrespondenceFilter::GetMTime()
{
  vtkMTimeType mTime = this->Superclass::GetMTime();
  vtkMTimeType time;

  // Would only make sense if it had inputs
  vtkTransform* transform = this->Transform;
  if (transform)
  {
    time = transform->GetMTime();
    mTime = time > mTime ? time : mTime;
  }
  return mTime;
}

int vtkPolyDataCorrespondenceFilter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkDebugMacro(<< __FUNCTION__);
  vtkPolyData* input0 = vtkPolyData::GetData(inputVector[0], 0);
  vtkPolyData* input1 = vtkPolyData::GetData(inputVector[1], 0);
  vtkPolyData* output0 = vtkPolyData::GetData(outputVector, 0);
  vtkPolyData* output1 = vtkPolyData::GetData(outputVector, 1);

  if (!input0 || !input1)
    return 0;

  if ((input0->GetPoints()->GetDataType() != input1->GetPoints()->GetDataType()))
  {
    vtkDebugMacro(<< "Input points must be of the same type");
    return 0;
  }

  // Initialize locators for the given search direction
  vtkPolyData* source = nullptr;
  vtkPolyData* target = nullptr;
  vtkImplicitPolyDataDistance2* locator;

  if (this->SearchDirection == SourceToTarget)
  {
    source = input0;
    locator = this->TargetLocator;
    target = input1;
  }
  else
  {
    source = input1;
    locator = this->SourceLocator;
    target = input0;
  }

  if (locator->GetInput())
  {
    if (locator->GetMTime() < source->GetMTime())
    {
      vtkDebugMacro(<< "Updating locator input");
      locator->SetInput(target);
    }
  }
  else
  {
    locator->SetInput(target);
  }

  int inputPrecision = source->GetPoints()->GetDataType();

  if (inputPrecision != VTK_FLOAT)
  {
    vtkDebugMacro(<< "We only support float input data");
    return 0;
  }

  vtkSmartPointer<vtkDataArray> sourceNormals = nullptr;

  // If we need source normals
  if (this->IncludeOriginNormals)
  {
    this->SourceNormals->SetInputData(source);

    // Unfortunately, it is only the extra points that respond to this
    this->SourceNormals->SetOutputPointsPrecision(this->OutputPrecision);
    this->SourceNormals->Update();
    sourceNormals = this->SourceNormals->GetOutput()->GetPointData()->GetNormals();

    if (sourceNormals->GetNumberOfTuples() != source->GetNumberOfPoints())
    {
      vtkDebugMacro(<< "Number of vertex normals does not match the number of points");
      return 0;
    }
  }

  vtkIdType step = 1;
  if (source->GetNumberOfPoints() > this->MaximumNumberOfLandmarks)
  {
    step = (source->GetNumberOfPoints() / this->MaximumNumberOfLandmarks);
  }

  vtkDebugMacro(<< "Executing functors");

  int outputPrecisionType = (this->OutputPrecision == vtkAlgorithm::DEFAULT_PRECISION
      ? source->GetPoints()->GetDataType()
      : (this->OutputPrecision == vtkAlgorithm::SINGLE_PRECISION ? VTK_FLOAT : VTK_DOUBLE));

  vtkNew<vtkPoints> points;
  vtkNew<vtkPoints> landmarks;
  vtkSmartPointer<vtkDataArray> landmarksNormals;
  vtkSmartPointer<vtkDataArray> originNormals = nullptr;

  if (outputPrecisionType == VTK_FLOAT)
  {
    landmarksNormals = vtkSmartPointer<vtkFloatArray>::New();
    originNormals = vtkSmartPointer<vtkFloatArray>::New();
    points->SetDataTypeToFloat();
    landmarks->SetDataTypeToFloat();
    if (this->IncludeOriginNormals)
    {
      NormalsFunctor<vtkFloatArray, vtkFloatArray, vtkFloatArray, vtkFloatArray> fun0(
        vtkFloatArray::SafeDownCast(source->GetPoints()->GetData()),
        vtkFloatArray::SafeDownCast(sourceNormals.GetPointer()),
        vtkFloatArray::SafeDownCast(points->GetData()),             // output
        vtkFloatArray::SafeDownCast(originNormals.GetPointer()),    // output
        vtkFloatArray::SafeDownCast(landmarks->GetData()),          // output
        vtkFloatArray::SafeDownCast(landmarksNormals.GetPointer()), // output
        step, locator, this);
      vtkSMPTools::For(0, source->GetNumberOfPoints(), fun0);
    }
    else
    {
      Functor<vtkFloatArray, vtkFloatArray, vtkFloatArray> fun(
        vtkFloatArray::SafeDownCast(source->GetPoints()->GetData()),
        vtkFloatArray::SafeDownCast(points->GetData()),             // output
        vtkFloatArray::SafeDownCast(landmarks->GetData()),          // output
        vtkFloatArray::SafeDownCast(landmarksNormals.GetPointer()), // output
        step, locator, this);
      vtkSMPTools::For(0, source->GetNumberOfPoints(), fun);
    }
  }
  else
  {
    landmarksNormals = vtkSmartPointer<vtkDoubleArray>::New();
    originNormals = vtkSmartPointer<vtkDoubleArray>::New();
    points->SetDataTypeToDouble();
    landmarks->SetDataTypeToDouble();
    Functor<vtkFloatArray, vtkDoubleArray, vtkDoubleArray> fun(
      vtkFloatArray::SafeDownCast(source->GetPoints()->GetData()),
      vtkDoubleArray::SafeDownCast(points->GetData()),
      vtkDoubleArray::SafeDownCast(landmarks->GetData()),
      vtkDoubleArray::SafeDownCast(landmarksNormals.GetPointer()), step, locator, this);
    vtkSMPTools::For(0, source->GetNumberOfPoints(), fun);
  }

  vtkDebugMacro(<< "Assigning output");
  landmarksNormals->SetName("Normals");
  output0->SetPoints(points);

  if (this->IncludeOriginNormals)
  {
    originNormals->SetName("Normals");
    output0->GetPointData()->AddArray(originNormals);
    output0->GetPointData()->SetActiveVectors("Normals");
  }
  output1->SetPoints(landmarks);
  output1->GetPointData()->AddArray(landmarksNormals);
  output1->GetPointData()->SetActiveVectors("Normals");
  return 1;
}
