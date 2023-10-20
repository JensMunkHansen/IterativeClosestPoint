#pragma once

#include "icp.hpp"

namespace icp
{
template <typename PointType, typename WeightType, int N>
class PointToPoint : public ICP<PointType, WeightType, N>
{
public:
  PointToPoint();
  ~PointToPoint();

  using typename ICP<PointType, WeightType, N>::AffineN;
  using typename ICP<PointType, WeightType, N>::MatrixNX;
  using typename ICP<PointType, WeightType, N>::VectorX;
  using typename ICP<PointType, WeightType, N>::VectorN;
  using typename ICP<PointType, WeightType, N>::MatrixNN;
  using typename ICP<PointType, WeightType, N>::MatrixXX;

  AffineN Update(MatrixNX& source, MatrixNX& target, const VectorX& weights);

  AffineN Update(
    MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights);

private:
  void ClosestRotation(const MatrixNX& M, MatrixNX& R);
};
} // namespace icp
