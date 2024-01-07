#pragma once

#include "icp.h"

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

  /*
   * Zhang, Zhengyou (1994). "Iterative point matching for registration of
   * free-form curves and surfaces". International Journal of Computer
   * Vision. 13 (12): 119–152. CiteSeerX 10.1.1.175.770.
   */
  AffineN Update(
    MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights);
  AffineN Update(MatrixNX& source, MatrixNX& target, const VectorX& weights);
  AffineN Update(MatrixNX& source, MatrixNX& target, const MatrixNX& sourceNormals,
    const MatrixNX& targetNormals, const VectorX& weights);

private:
};
} // namespace icp
