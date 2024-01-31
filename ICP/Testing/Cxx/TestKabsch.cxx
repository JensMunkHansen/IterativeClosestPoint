#include <icp/kabsch.h>
#include <vtkMultiMeshAlignment.h>
#include <vtkNew.h>
#include <vtkSetGet.h>

int TestKabsch(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  vtkNew<vtkMultiMeshAlignment> alignment;
  alignment->TestKabsch();
  return 0;
}
