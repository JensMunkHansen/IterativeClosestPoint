#include <icp/plane2plane.hpp>

namespace icp
{
template <typename PointType, typename WeightType, int N>
PlaneToPlane<PointType, WeightType, N>::PlaneToPlane()
{
  // Point-to-plane only works in 3D
  static_assert(N == 3, "Point-to-plane only works in 3D!");
}

template <typename PointType, typename WeightType, int N>
PlaneToPlane<PointType, WeightType, N>::~PlaneToPlane()
{
}

template <typename PointType, typename WeightType, int N>
typename PlaneToPlane<PointType, WeightType, N>::AffineN
PlaneToPlane<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const VectorX& weights)
{
  AffineN transform;
  return transform;
}

template <typename PointType, typename WeightType, int N>
typename PlaneToPlane<PointType, WeightType, N>::AffineN
PlaneToPlane<PointType, WeightType, N>::Update(
  MatrixNX& source, MatrixNX& target, const MatrixNX& targetNormals, const VectorX& weights)
{
  AffineN transform;
  return transform;
}

template <typename PointType, typename WeightType, int N>
typename PlaneToPlane<PointType, WeightType, N>::AffineN
PlaneToPlane<PointType, WeightType, N>::Update(
  MatrixNX& P, MatrixNX& Q, const MatrixNX& NP, const MatrixNX& NQ, const VectorX& weights)
{
  // normalize point sets
  VectorN Pmean = P.rowwise().mean();
  VectorN Qmean = Q.rowwise().mean();

  MatrixNX Pbar = P.colwise() - Pmean;
  MatrixNX Qbar = Q.colwise() - Qmean;

  // sum of normals
  MatrixNX Norm = NP + NQ;

  // compute A and b of linear system
  int num_points = P.cols();

  Matrix66 LHS = Matrix66::Zero();
  Vector6 RHS = Vector6::Zero();

  for (int i = 0; i < num_points; ++i)
  {
    VectorN n_i = Norm.col(i);
    VectorN p_i = Pbar.col(i);
    VectorN q_i = Qbar.col(i);
    typename MatrixXX::Scalar b_i = (p_i - q_i).dot(n_i);

    Vector6 x_i = Vector6::Zero();
    x_i << (p_i + q_i).cross(n_i), n_i;
    LHS += x_i * x_i.transpose();
    RHS += b_i * x_i;
  }
  LHS = LHS.template selfadjointView<vtkeigen::Upper>();
  vtkeigen::LDLT<Matrix66> ldlt(LHS);
  RHS = ldlt.solve(-RHS).transpose();

  VectorN a_tilda = RHS.template head<3>();
  VectorN t_tilda = RHS.template tail<3>();

  // compute intermediate rotation
  typename MatrixXX::Scalar theta = atan(a_tilda.norm()); // rotation angle
  VectorN a = a_tilda.normalized();                       // normalized rotation axis
  MatrixNN W;
  W << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;
  MatrixNN intermediate_R = MatrixNN::Identity() + sin(theta) * W + (1 - cos(theta)) * (W * W);

  // compose translations and rotations
  VectorN t1 = -Pmean.transpose();
  VectorN t2 = cos(theta) * t_tilda;
  VectorN t3 = Qmean.transpose();
  MatrixNN R = intermediate_R * intermediate_R;
  vtkeigen::Matrix<PointType, 1, 3> t =
    (intermediate_R * intermediate_R * t1) + (intermediate_R * t2) + t3;

  AffineN transform;
  transform.linear() = R;
  transform.translation() = t;
  return transform;
}
#ifdef ICP_DOUBLE_SUPPORT
template class PlaneToPlane<double, double, 3U>;
#endif
template class PlaneToPlane<float, float, 3U>;

}
