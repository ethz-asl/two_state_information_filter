#ifndef TSIF_COMMON_HPP_
#define TSIF_COMMON_HPP_

#include <array>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>

#if TSIF_VERBOSE > 0
#define TSIF_LOG(msg) std::cout << msg << std::endl
#define TSIF_LOGIF(con,msg) if(con) TSIF_LOG(msg)
#define TSIF_LOGW(msg) std::cout << "\033[33m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGWIF(con,msg) if(con) TSIF_LOGW(msg)
#define TSIF_LOGE(msg) std::cout << "\033[31m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGEIF(con,msg) if(con) TSIF_LOGE(msg)
#else
#define TSIF_LOG(msg)
#define TSIF_LOGIF(con,msg)
#define TSIF_LOGW(msg)
#define TSIF_LOGWIF(con,msg)
#define TSIF_LOGE(msg)
#define TSIF_LOGEIF(con,msg)
#endif

namespace tsif{

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

typedef Eigen::Quaterniond Quat;
static constexpr double Sinc(const double x){
  return fabs(x) < 1e-8 ? 1 : sin(x)/x;
}
static Quat Exp(const Vec3& v){
  const double ha = 0.5*v.norm();
  const double re = cos(ha);
  const Vec3 im = 0.5*Sinc(ha)*v;
  return Quat(re,im(0),im(1),im(2));
}
static Vec3 Log(const Quat& q){
  const double re = q.w();
  const Vec3 im(q.x(),q.y(),q.z());
  const double sha = im.norm();
  return sha < 1e-8 ? ((std::signbit(re)!=0)*-2+1)*2*im : 2*atan2(sha,re)/sha*im;
}
static Quat Boxplus(const Quat& q, const Vec3& v){
  return Exp(v)*q;
}
static Vec3 Boxminus(const Quat& q, const Quat& p){
  return Log(q*p.inverse());
}
static Mat3 SSM(const Vec3& vec){
  Mat3 mat;
  mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
  return mat;
}
static Mat3 RotMat(const Vec3& vec){
  const double a = vec.norm();
  if(a < 1e-8){
    return Mat3::Identity() + SSM(vec);
  } else{
    const Vec3 axis = vec.normalized();
    return Mat3::Identity() + sin(a)*SSM(axis) + (1-cos(a))*SSM(axis)*SSM(axis);
  }
}
static Mat3 GammaMat(const Vec3& vec){
  const double a = vec.norm();
  if(a < 1e-8){
    return Mat3::Identity() + 0.5 * SSM(vec);
  } else{
    const Vec3 axis = vec.normalized();
    return Mat3::Identity() + (1-cos(a))/a*SSM(axis) + (a-sin(a))/a*SSM(axis)*SSM(axis);
  }
}
static Mat3 FromTwoVectorsJac(const Vec3& a, const Vec3& b){
  const Vec3 cross = a.cross(b);
  const double crossNorm = cross.norm();
  Vec3 crossNormalized = cross/crossNorm;
  Mat3 crossNormalizedSqew = SSM(crossNormalized);
  const double c = a.dot(b);
  const double angle = std::acos(c);
  if(crossNorm<1e-6){
    if(c>0){
      return -SSM(b);
    } else {
      TSIF_LOGW("Warning: instable FromTwoVectorsJac!");
      return Mat3::Zero();
    }
  } else {
    return -1/crossNorm*(crossNormalized*b.transpose()-(crossNormalizedSqew*crossNormalizedSqew*SSM(b)*angle));
  }
}

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double sec){
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration){
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}
static std::string Print(TimePoint t){
  std::ostringstream out;
  out.precision(15);
  out << ((double)t.time_since_epoch().count()*Duration::period::num)/(Duration::period::den);
  return out.str();
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
  static NormalRandomNumberGenerator& Instance(){
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

} // namespace tsif

#endif /* TSIF_COMMON_HPP_ */
