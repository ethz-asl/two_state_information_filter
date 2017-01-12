#ifndef GIF_COMMON_HPP_
#define GIF_COMMON_HPP_

#include <array>
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "kindr/Core"
#include <chrono>
#include <glog/logging.h>

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
using MatRef = typename std::conditional<(N==1 & M>1),
                                         Eigen::Ref<Mat<N,M>,0,Eigen::InnerStride<>>,
                                         Eigen::Ref<Mat<N,M>>>::type;
using MatRef3 = MatRef<3>;
using MatRefX = MatRef<>;
template<int N = -1, int M = N>
using MatCRef = Eigen::Ref<const Mat<N,M>>;
using MatCRef3 = MatCRef<3>;
using MatCRefX = MatCRef<>;

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double sec) {
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}
static std::string Print(TimePoint t){
  std::ostringstream out;
  out.precision(15);
  out << ((double)t.time_since_epoch().count()*Duration::period::num)/(Duration::period::den);
  return out.str();
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

/*! \brief Normal Random Number Generator
 *         Singleton class for generating normal random numbers (N(0,1)). Allows setting of seed.
 */
class NormalRandomNumberGenerator
{
 public:
  void SetSeed(int s){
    generator_.seed(s);
  }
  double Get(){
    return distribution_(generator_);
  }
  template<int N>
  Vec<N> GetVec(){
    Vec<N> n;
    for(int i=0;i<N;i++){
      n(i) = Get();
    }
    return n;
  }
  static NormalRandomNumberGenerator& Instance()
  {
    static NormalRandomNumberGenerator instance;
    return instance;
  }
  std::default_random_engine& GetGenerator(){
    return generator_;
  }
 protected:
  std::default_random_engine generator_;
  std::normal_distribution<double> distribution_;
  NormalRandomNumberGenerator(): generator_(0), distribution_(0.0,1.0){
  }
};

}

#endif /* GIF_COMMON_HPP_ */
