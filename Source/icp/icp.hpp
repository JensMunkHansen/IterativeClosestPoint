#pragma once

#include <icp/config.h>
#include "vtk_eigen.h"
#include VTK_EIGEN(Dense)

namespace icp
{
template <typename PointType, typename WeightType, int N>
class ICP
{
public:
  // Input/output types
  using AffineN = vtkeigen::Transform<PointType, N, vtkeigen::Affine>;
  using MatrixNX = vtkeigen::Matrix<PointType, N, vtkeigen::Dynamic>;
  using MatrixXX = vtkeigen::Matrix<PointType, vtkeigen::Dynamic, vtkeigen::Dynamic>;
  using VectorX = vtkeigen::Matrix<PointType, vtkeigen::Dynamic, 1>;
  using VectorN = vtkeigen::Matrix<PointType, N, 1>;
  using MatrixNN = vtkeigen::Matrix<PointType, N, N>;

  // Special for point-to-plane in 3D
  using Matrix66 = vtkeigen::Matrix<PointType, 6, 6>;
  using Vector6 = vtkeigen::Matrix<PointType, 6, 1>;
  using Vector3 = vtkeigen::Matrix<PointType, 3, 1>;
  using Block33 = vtkeigen::Block<Matrix66, 3, 3>;

  virtual AffineN Update(
    MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights) = 0;
  virtual AffineN Update(MatrixNX& source, MatrixNX& target, const VectorX& weights) = 0;

  virtual ~ICP() = default;

protected:
};
}  // namespace icp
