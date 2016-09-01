#ifndef GIF_COMMON_HPP_
#define GIF_COMMON_HPP_

#include <array>
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "kindr/Core"
#include <chrono>

namespace GIF {

typedef kindr::RotationQuaternionPD Quat;
typedef kindr::RotationMatrixPD RotMat;

template<int N = -1>
using Vec = Eigen::Matrix<double,N,1>;
using Vec3 = Vec<3>;
using VecX = Vec<>;
template<int N = -1>
using VecRef = Eigen::Ref<Vec<N>>;
using VecRef3 = VecRef<3>;
using VecRefX = VecRef<>;
template<int N = -1>
using VecCRef = Eigen::Ref<const Vec<N>>;
using VecCRef3 = VecCRef<3>;
using VecCRefX = VecCRef<>;

template<int N = -1, int M = N>
using Mat = Eigen::Matrix<double,N,M>;
using Mat3 = Mat<3>;
using MatX = Mat<>;
template<int N = -1, int M = N>
using MatRef = Eigen::Ref<Mat<N,M>>;
template<int N = -1, int M = N>
using MatCRef = Eigen::Ref<const Mat<N,M>>;

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double sec) {
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration)
      .count();
}
inline Mat3 gSM(const Vec3& vec) {
  return kindr::getSkewMatrixFromVector(vec);
}

static void enforceSymmetry(MatX& mat) {
  mat = 0.5 * (mat + mat.transpose()).eval();
}

inline Mat3 GammaMat(const Vec3& a) {
  return kindr::getJacobianOfExponentialMap(a);
}

}

#endif /* GIF_COMMON_HPP_ */
