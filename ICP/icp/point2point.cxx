#include "point2point.hpp"
#include <iostream>
namespace icp
{
template <typename PointType, typename WeightType, int N>
PointToPoint<PointType, WeightType, N>::PointToPoint()
{
}

template <typename PointType, typename WeightType, int N>
PointToPoint<PointType, WeightType, N>::~PointToPoint()
{
}

template <typename PointType, typename WeightType, int N>
typename PointToPoint<PointType, WeightType, N>::AffineN
PointToPoint<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const VectorX& weights)
{
  /// Normalize weight vector
  VectorX w_normalized = weights / weights.sum();
  /// De-mean
  VectorX sourceMean(N), targetMean(N);

  for (int i = 0; i < N; ++i)
  {
    sourceMean(i) = (source.row(i).array() * w_normalized.transpose().array()).sum();
    targetMean(i) = (target.row(i).array() * w_normalized.transpose().array()).sum();
  }
  source.colwise() -= sourceMean;
  target.colwise() -= targetMean;

  /// Compute transformation
  AffineN transformation;

  MatrixNN sigma = source * w_normalized.asDiagonal() * target.transpose();
  Eigen::JacobiSVD<MatrixNN> svd(sigma, vtkeigen::ComputeFullU | vtkeigen::ComputeFullV);
  if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0.0)
  {
    // Strict rotations - no reflections
    VectorN S = VectorN::Ones(N);
    S(N - 1) = -1.0;
    transformation.linear() = svd.matrixV() * S.asDiagonal() * svd.matrixU().transpose();
  }
  else
  {
    transformation.linear() = svd.matrixV() * svd.matrixU().transpose();
  }

  transformation.translation() = targetMean - transformation.linear() * sourceMean;

  /// Re-apply mean (no need, really)
  source.colwise() += sourceMean;
  target.colwise() += targetMean;

  /// Return transformation
  return transformation;
}

template <typename PointType, typename WeightType, int N>
typename PointToPoint<PointType, WeightType, N>::AffineN
PointToPoint<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights)
{
  return this->Update(source, target, weights);
}

#ifdef ICP_DOUBLE_SUPPORT
template class PointToPoint<double, double, 3U>;
template class PointToPoint<double, double, 2U>;
#endif
template class PointToPoint<float, float, 3U>;
template class PointToPoint<float, float, 2U>;
}
