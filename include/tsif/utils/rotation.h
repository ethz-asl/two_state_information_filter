#ifndef TSIF_ROTATION_HPP_
#define TSIF_ROTATION_HPP_

#include "tsif/utils/logging.h"
#include "tsif/utils/typedefs.h"

namespace tsif{

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
  return sha < 1e-8 ? ((std::signbit(re)!=0)*-2+1)*2*im : 2*atan(sha/re)/sha*im;
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

} // namespace tsif

#endif /* TSIF_ROTATION_HPP_ */
