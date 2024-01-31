#pragma once

#include "vtk_eigen.h"
#include <icp/config.h>
#include VTK_EIGEN(Dense)
#include <vector>

namespace icp
{
template <typename PointType>
class Kabsch
{
public:
  Kabsch();
  ~Kabsch();
  using Vector3 = vtkeigen::Matrix<PointType, 3, 1>;
  using Matrix3 = vtkeigen::Matrix<PointType, 3, 3>;
  using Affine3 = vtkeigen::Transform<PointType, 3, vtkeigen::Affine>;
  using MatrixX = vtkeigen::Matrix<PointType, vtkeigen::Dynamic, vtkeigen::Dynamic>;
  bool Update(const std::vector<Vector3>& target, const std::vector<Vector3>& source, Affine3& tx,
    PointType& rms);
};

} // namespace icp
