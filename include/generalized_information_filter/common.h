#ifndef GIF_COMMON_HPP_
#define GIF_COMMON_HPP_

#include <array>
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "kindr/Core"
#include <chrono>

namespace GIF {

typedef kindr::RotationQuaternionPD QPD;
typedef kindr::RotationMatrixPD MPD;
typedef Eigen::Vector3d V3D;
typedef Eigen::Matrix3d M3D;
typedef Eigen::VectorXd VXD;
typedef Eigen::MatrixXd MXD;

template<int N = -1>
using Vec = Eigen::Matrix<double,N,1>;
template<int N = -1>
using VecR = Eigen::Ref<Vec<N>>;
template<int N = -1>
using VecRC = Eigen::Ref<const Vec<N>>;

template<int N = -1, int M = N>
using Mat = Eigen::Matrix<double,N,M>;
template<int N = -1, int M = N>
using MatR = Eigen::Ref<Mat<N,M>>;
template<int N = -1, int M = N>
using MatRC = Eigen::Ref<const Mat<N,M>>;

template<typename T>
using SP = std::shared_ptr<T>;

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double& sec) {
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration)
      .count();
}
inline M3D gSM(const V3D& vec) {
  return kindr::getSkewMatrixFromVector(vec);
}

static void enforceSymmetry(MXD& mat) {
  mat = 0.5 * (mat + mat.transpose()).eval();
}

inline M3D Lmat(const V3D& a) {
  return kindr::getJacobianOfExponentialMap(a);
}

}

#endif /* GIF_COMMON_HPP_ */
