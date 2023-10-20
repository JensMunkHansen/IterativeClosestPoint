#include "vtkPolyDataCorrespondenceFilter.h"

#include <vtkParametricFunctionSource.h>
#include <vtkParametricRandomHills.h>
#include <vtkPointData.h>
#include <vtkSMPTools.h>
#include <vtkTimerLog.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

namespace
{
vtkSmartPointer<vtkPolyDataAlgorithm> CreateRandomHills(int nHills = 20)
{
  vtkNew<vtkParametricRandomHills> sourceFcn;
  sourceFcn->AllowRandomGenerationOn();
  sourceFcn->SetRandomSeed(1);
  sourceFcn->SetNumberOfHills(nHills);
  vtkNew<vtkParametricFunctionSource> source;

  source->SetParametricFunction(vtkParametricFunction::SafeDownCast(sourceFcn));
  source->SetUResolution(51);
  source->SetVResolution(51);
  source->SetWResolution(51);
  source->Update();
  return source;
}
}

int TestPolyDataCorrespondenceFilter(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  vtkSMPTools::SetBackend("TBB");
  // vtkSMPTools::SetBackend("Sequential");

  vtkSmartPointer<vtkPolyDataAlgorithm> sourceFilter = CreateRandomHills();
  vtkPolyData* source = vtkPolyData::SafeDownCast(sourceFilter->GetOutput());

  vtkNew<vtkTransform> transform;
  transform->Translate(0.0, 0.0, -3.0);
  vtkNew<vtkTransformPolyDataFilter> transformPoly;
  transformPoly->SetTransform(transform);
  transformPoly->SetInputConnection(sourceFilter->GetOutputPort());
  transformPoly->Update();
  vtkPolyData* target = transformPoly->GetOutput();

  std::cout << "Number of target points: " << target->GetNumberOfPoints() << std::endl;

  vtkNew<vtkPolyDataCorrespondenceFilter> correspondences;
  correspondences->SetMaximumNumberOfLandmarks(20000);
  correspondences->SetInputData(0, source);
  correspondences->SetInputData(1, target);
  // correspondences->SetOutputPrecision(vtkAlgorithm::DOUBLE_PRECISION);
  correspondences->SetOutputPrecision(vtkAlgorithm::SINGLE_PRECISION);
  //  correspondences->SetIncludeSourceNormals(true);
  vtkNew<vtkTimerLog> tl;
  tl->StartTimer();
  for (int i = 0; i < 100; i++)
  {
    correspondences->Update();
    correspondences->Modified();
  }
  tl->StopTimer();
  double elapsed = tl->GetElapsedTime() / 100;
  std::cout << "Elapsed: " << elapsed << std::endl;

  std::cout << "Number of origins: " << correspondences->GetOutput(0)->GetNumberOfPoints()
            << std::endl;

  std::cout << "Number of landmarks: " << correspondences->GetOutput(1)->GetNumberOfPoints()
            << std::endl;
  std::cout << "Data type (origins): " << correspondences->GetOutput(0)->GetPoints()->GetDataType()
            << std::endl;
  std::cout << "Data type (landmarks): "
            << correspondences->GetOutput(1)->GetPoints()->GetDataType() << std::endl;
  std::cout << "Data type (landmark normals): "
            << correspondences->GetOutput(1)->GetPointData()->GetArray("Normals")->GetDataType()
            << std::endl;
  vtkIdType nOriginNormals = 0;
  vtkDataArray* originNormals = correspondences->GetOutput(0)->GetPointData()->GetArray("Normals");
  if (originNormals)
  {
    nOriginNormals = originNormals->GetNumberOfTuples();
  }
  std::cout << "Number of origin normals: " << nOriginNormals << std::endl;
  return 0;
}
