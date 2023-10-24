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

  /*
    Besl, Paul J.; N.D. McKay (1992). "A Method for Registration of
    3-D Shapes". IEEE Transactions on Pattern Analysis and Machine
    Intelligence. 14 (2): 239–256. doi:10.1109/34.121791
  */
  AffineN Update(MatrixNX& source, MatrixNX& target, const VectorX& weights);

  AffineN Update(
    MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights);

  AffineN Update(MatrixNX& source, MatrixNX& target, const MatrixNX& sourceNormals,
    const MatrixNX& targetNormals, const VectorX& weights);

private:
};
} // namespace icp
