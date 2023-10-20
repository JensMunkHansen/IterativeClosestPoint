#pragma once

#include "icp.hpp"

namespace icp
{
template <typename PointType, typename WeightType, int N>
class PointToPlane : public ICP<PointType, WeightType, N>
{
public:
  PointToPlane();
  ~PointToPlane();

  using typename ICP<PointType, WeightType, N>::AffineN;
  using typename ICP<PointType, WeightType, N>::MatrixNX;
  using typename ICP<PointType, WeightType, N>::VectorX;
  using typename ICP<PointType, WeightType, N>::VectorN;
  using typename ICP<PointType, WeightType, N>::MatrixNN;
  using typename ICP<PointType, WeightType, N>::MatrixXX;

  using typename ICP<PointType, WeightType, N>::Vector3;
  using typename ICP<PointType, WeightType, N>::Vector6;
  using typename ICP<PointType, WeightType, N>::Matrix66;
  using typename ICP<PointType, WeightType, N>::Block33;

#if 0
  AffineN Symmetric(
    MatrixNX& P, MatrixNX& Q, const MatrixNX& NP, const MatrixNX& NQ, const VectorX& W);
#endif

  AffineN Update(
    MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights);
  AffineN Update(MatrixNX& source, MatrixNX& target, const VectorX& weights);

private:
};
} // namespace icp
