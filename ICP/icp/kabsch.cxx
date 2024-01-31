#include <icp/kabsch.h>

namespace icp
{

/*
 * Estimate _tx such that:  _tgt = _tx * _src
 * Requires at least four pairs of points.
 * See:
 *   https://en.wikipedia.org/wiki/Kabsch_algorithm
 *   https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
 *   Closed-form solution of absolute orientation using unit quaternions - Horn 1987
 */

template <typename PointType>
Kabsch<PointType>::Kabsch()
{
}
template <typename PointType>
Kabsch<PointType>::~Kabsch()
{
}

template <typename PointType>
bool Kabsch<PointType>::Update(const std::vector<Vector3>& _tgt, const std::vector<Vector3>& _src,
  Affine3& _tx, PointType& _residual_rms)
{
  if (_tgt.size() != _src.size())
  {
    printf("Kabsch: _tgt.size() != _src.size() (%i and %i respectively). Check input data\n",
      _tgt.size(), _src.size());
    return false;
  }
  if (_tgt.size() < 4)
  {
    printf(
      "Kabsch: need at least four point pairs (more is better). Only got %i points\n", _tgt.size());
    return false;
  }

  size_t n = _tgt.size();
  MatrixX src(3, n), tgt(3, n);
  for (int i = 0; i < n; i++)
  {
    tgt.col(i) = _tgt[i];
    src.col(i) = _src[i];
  }

  Vector3 centroid_src(0, 0, 0), centroid_tgt(0, 0, 0);
  for (int c = 0; c < n; c++)
  {
    centroid_src += src.col(c);
    centroid_tgt += tgt.col(c);
  }
  centroid_src /= n;
  centroid_tgt /= n;
  for (int c = 0; c < n; c++)
  {
    src.col(c) -= centroid_src;
    tgt.col(c) -= centroid_tgt;
  }

  PointType scale = 1.0;
  PointType mean_dist_src = 0, mean_dist_tgt = 0;
  {
    for (int c = 0; c < n; c++)
    {
      mean_dist_src += src.col(c).norm();
      mean_dist_tgt += tgt.col(c).norm();
    }
    printf("Mean dist src=%.12f\n", mean_dist_src);
    printf("Mean dist target=%.12f\n", mean_dist_tgt);

    if (std::fabs(mean_dist_tgt - mean_dist_src) < 10 * std::numeric_limits<PointType>::epsilon())
    {
      scale = PointType(1.0);
    }
    else
    {
      scale = mean_dist_tgt / mean_dist_src;
    }
    printf("Kabsch: scale=%.12f\n", scale);
    centroid_tgt /= scale;
  }

  /*
   * Kabsch estimates the rotation matrix which minimizes the RMSD.
   * The error is a function of both the angle and the length of the
   * vector - that is the distance from the centroid.
   *
   * This selective normalization avoids biasing the result towards
   * points which are further than the average distance from the centroid.
   */

  for (int c = 0; c < n; c++)
  {
    if (src.col(c).norm() > mean_dist_src)
      src.col(c).normalize();
    else
      src.col(c) /= std::max<PointType>(mean_dist_src, std::numeric_limits<PointType>::epsilon());

    if (tgt.col(c).norm() > mean_dist_tgt)
      tgt.col(c).normalize();
    else
      tgt.col(c) /= std::max<PointType>(mean_dist_tgt, std::numeric_limits<PointType>::epsilon());
  }

  Matrix3 rotation = Matrix3::Identity();
  {
    MatrixX cov = src * tgt.transpose();
    Eigen::JacobiSVD<MatrixX> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    PointType d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
      d = 1.0;
    else
      d = -1.0;
    Matrix3 I = Matrix3::Identity();
    I(2, 2) = d;
    rotation = svd.matrixV() * I * svd.matrixU().transpose();
  }

  _tx.setIdentity();
  _tx.linear() = scale * rotation;
  _tx.translation() = scale * (centroid_tgt - rotation * centroid_src);

  _residual_rms = 0;
  for (int c = 0; c < n; c++)
  {
    auto err = (_tgt[c] - _tx * _src[c]).norm();
    _residual_rms += err * err;
  }
  _residual_rms /= n;
  _residual_rms = sqrt(_residual_rms);

  printf("Kabsch: RMSD=%.32f\n", _residual_rms);

  return true;
}
#ifdef ICP_DOUBLE_SUPPORT
template class Kabsch<double>;
#endif
template class Kabsch<float>;
}
