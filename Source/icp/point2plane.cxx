#include "point2plane.hpp"
#include <iostream>

namespace icp
{
template <typename PointType, typename WeightType, int N>
PointToPlane<PointType, WeightType, N>::PointToPlane()
{
  // Point-to-plane only works in 3D
  static_assert(N == 3, "Point-to-plane only works in 3D!");
}

template <typename PointType, typename WeightType, int N>
PointToPlane<PointType, WeightType, N>::~PointToPlane()
{
}

template <typename PointType, typename WeightType, int N>
typename PointToPlane<PointType, WeightType, N>::AffineN
PointToPlane<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const VectorX& weights)
{
  AffineN transform;
  return transform;
}

template <typename PointType, typename WeightType, int N>
typename PointToPlane<PointType, WeightType, N>::AffineN
PointToPlane<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights)
{
  /// Normalize weight vector
  VectorX weightNormalized = weights / weights.sum();
  // De-mean
  Vector3 sourceMean;
  for (int i = 0; i < N; i++)
    sourceMean(i) = (source.row(i).array() * weightNormalized.transpose().array()).sum();
  source.colwise() -= sourceMean;
  target.colwise() -= sourceMean;

  // Prepare LHS and RHS for LDLT:
  Matrix66 LHS = Matrix66::Zero();
  Vector6 RHS = Vector6::Zero();
  Block33 TL = LHS.template topLeftCorner<N, N>();
  Block33 TR = LHS.template topRightCorner<N, N>();
  Block33 BR = LHS.template bottomRightCorner<N, N>();

  // Why is this okay? <dynamic, N> = <dynamic, dynamic>
  MatrixNX C = MatrixXX::Zero(N, source.cols());

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < source.cols(); i++)
    {
      C.col(i) = source.col(i).cross(targetNormals.col(i));
    }
#pragma omp sections nowait
    {
#pragma omp section
      for (int i = 0; i < source.cols(); i++)
        TL.template selfadjointView<vtkeigen::Upper>().rankUpdate(C.col(i), weights(i));
#pragma omp section
      for (int i = 0; i < source.cols(); i++)
        TR += (C.col(i) * targetNormals.col(i).transpose()) * weights(i);
#pragma omp section
      for (int i = 0; i < source.cols(); i++)
        BR.template selfadjointView<vtkeigen::Upper>().rankUpdate(targetNormals.col(i), weights(i));
#pragma omp section
      for (int i = 0; i < C.cols(); i++)
      {
        PointType distToPlane =
          -((source.col(i) - target.col(i)).dot(targetNormals.col(i))) * weights(i);
        RHS.template head<3>() += C.col(i) * distToPlane;
        RHS.template tail<3>() += targetNormals.col(i) * distToPlane;
      }
    }
  }

  LHS = LHS.template selfadjointView<vtkeigen::Upper>();

  /// Compute transformation
  AffineN transformation;

  vtkeigen::LDLT<Matrix66> ldlt(LHS);
  RHS = ldlt.solve(RHS).transpose();

  // Point-to-plane does not make sense in 2D
  transformation = vtkeigen::AngleAxis<PointType>(RHS(0), VectorN::UnitX()) *
    vtkeigen::AngleAxis<PointType>(RHS(1), VectorN::UnitY()) *
    vtkeigen::AngleAxis<PointType>(RHS(2), VectorN::UnitZ());

  transformation.translation() = RHS.template tail<N>();

  /// Re-apply mean
  source.colwise() += sourceMean;
  target.colwise() += sourceMean;
  transformation.translation() += sourceMean - transformation.linear() * sourceMean;

  /// Return transformation
  return transformation;
}
#ifdef ICP_DOUBLE_SUPPORT
  template class PointToPlane<double, double, 3U>;
#endif
template class PointToPlane<float, float, 3U>;
}
