#include "tsif/utils/camera.h"

namespace tsif{

  Camera::Camera(){
    k1_ = -0.28340811; k2_ = 0.07395907; k3_ = 0.0; k4_ = 0.0; k5_ = 0.0; k6_ = 0.0;
    p1_ = 0.00019359; p2_ = 1.76187114e-05; s1_ = 0.0; s2_ = 0.0; s3_ = 0.0; s4_ = 0.0;
    K_.setIdentity();
    K_(0,0) = 458.654;
    K_(0,2) = 367.215;
    K_(1,1) = 457.296;
    K_(1,2) = 248.375;
    type_ = RADTAN;
  };

  Camera::~Camera(){};

  void Camera::LoadCameraMatrix(const std::string& filename){
  }

  void Camera::LoadRadtan(const std::string& filename){
  }

  void Camera::LoadEquidist(const std::string& filename){
  }

  void Camera::Load(const std::string& filename){
  }

  void Camera::DistortRadtan(const Eigen::Vector2d& in, Eigen::Vector2d& out) const{
    const double x2 = in(0) * in(0);
    const double y2 = in(1) * in(1);
    const double xy = in(0) * in(1);
    const double r2 = x2 + y2;
    const double kr = (1 + ((k3_ * r2 + k2_) * r2 + k1_) * r2);
    out(0) = in(0) * kr + p1_ * 2 * xy + p2_ * (r2 + 2 * x2);
    out(1) = in(1) * kr + p1_ * (r2 + 2 * y2) + p2_ * 2 * xy;
  }
  void Camera::DistortRadtan(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const{
    const double x2 = in(0) * in(0);
    const double y2 = in(1) * in(1);
    const double xy = in(0) * in(1);
    const double r2 = x2 + y2;
    const double kr = (1 + ((k3_ * r2 + k2_) * r2 + k1_) * r2);
    out(0) = in(0) * kr + p1_ * 2 * xy + p2_ * (r2 + 2 * x2);
    out(1) = in(1) * kr + p1_ * (r2 + 2 * y2) + p2_ * 2 * xy;
    J(0,0) = kr + 2.0 * k1_ * x2 + 4.0 * k2_ * x2 * r2 + 6.0 * k3_ * x2 * r2 * r2 + 2.0 * p1_ * in(1) + 6.0 * p2_ * in(0);
    J(0,1) = 2.0 * k1_ * xy + 4.0 * k2_ * xy * r2 + 6.0 * k3_ * xy * r2 * r2 + 2 * p1_ * in(0) + 2 * p2_ * in(1);
    J(1,0) = J(0,1);
    J(1,1) = kr + 2.0 * k1_ * y2 + 4.0 * k2_ * y2 * r2 + 6.0 * k3_ * y2 * r2 * r2 + 6.0 * p1_ * in(1) + 2.0 * p2_ * in(0);
  }

  void Camera::DistortEquidist(const Eigen::Vector2d& in, Eigen::Vector2d& out) const{
    const double x2 = in(0) * in(0);
    const double y2 = in(1) * in(1);
    const double r = std::sqrt(x2 + y2); // 1/r*x

    if(r < 1e-8){
      out(0) = in(0);
      out(1) = in(1);
      return;
    }

    const double th = atan(r); // 1/(r^2 + 1)
    const double th2 = th*th;
    const double th4 = th2*th2;
    const double th6 = th2*th4;
    const double th8 = th2*th6;
    const double thd = th * (1.0 + k1_ * th2 + k2_ * th4 + k3_ * th6 + k4_ * th8);
    const double s = thd/r;

    out(0) = in(0) * s;
    out(1) = in(1) * s;
  }
  void Camera::DistortEquidist(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const{
    const double x2 = in(0) * in(0);
    const double y2 = in(1) * in(1);
    const double r = std::sqrt(x2 + y2);

    if(r < 1e-8){
      out(0) = in(0);
      out(1) = in(1);
      J.setIdentity();
      return;
    }

    const double r_x = 1/r*in(0);
    const double r_y = 1/r*in(1);

    const double th = atan(r); // 1/(r^2 + 1)
    const double th_r = 1/(r*r+1);
    const double th2 = th*th;
    const double th4 = th2*th2;
    const double th6 = th2*th4;
    const double th8 = th2*th6;
    const double thd = th * (1.0 + k1_ * th2 + k2_ * th4 + k3_ * th6 + k4_ * th8);
    const double thd_th = 1.0 + 3 * k1_ * th2 + 5* k2_ * th4 + 7 * k3_ * th6 + 9 * k4_ * th8;
    const double s = thd/r;
    const double s_r = thd_th*th_r/r - thd/(r*r);

    out(0) = in(0) * s;
    out(1) = in(1) * s;

    J(0,0) = s + in(0)*s_r*r_x;
    J(0,1) = in(0)*s_r*r_y;
    J(1,0) = in(1)*s_r*r_x;
    J(1,1) = s + in(1)*s_r*r_y;
  }

  void Camera::Distort(const Eigen::Vector2d& in, Eigen::Vector2d& out) const{
    switch(type_){
      case RADTAN:
        DistortRadtan(in,out);
        break;
      case EQUIDIST:
        DistortEquidist(in,out);
        break;
      default:
        DistortRadtan(in,out);
        break;
    }
  }
  void Camera::Distort(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const{
    switch(type_){
      case RADTAN:
        DistortRadtan(in,out,J);
        break;
      case EQUIDIST:
        DistortEquidist(in,out,J);
        break;
      default:
        DistortRadtan(in,out,J);
        break;
    }
  }

  bool Camera::BearingToPixel(const Eigen::Vector3d& vec, Eigen::Vector2d& c) const{
    // Project
    if(vec(2)<=1e-8) return false;
    const Eigen::Vector2d undistorted = Eigen::Vector2d(vec(0)/vec(2),vec(1)/vec(2));

    // Distort
    Eigen::Vector2d distorted;
    Distort(undistorted,distorted);

    // Shift origin and scale
    c(0) = K_(0, 0)*distorted(0) + K_(0, 2);
    c(1) = K_(1, 1)*distorted(1) + K_(1, 2);
    return true;
  }
  bool Camera::BearingToPixel(const Eigen::Vector3d& vec, Eigen::Vector2d& c, Eigen::Matrix<double,2,3>& J) const{
    // Project
    if(vec(2)<=1e-8) return false;
    const Eigen::Vector2d undistorted = Eigen::Vector2d(vec(0)/vec(2),vec(1)/vec(2));
    Eigen::Matrix<double,2,3> J1; J1.setZero();
    J1(0,0) = 1.0/vec(2);
    J1(0,2) = -vec(0)/pow(vec(2),2);
    J1(1,1) = 1.0/vec(2);
    J1(1,2) = -vec(1)/pow(vec(2),2);

    // Distort
    Eigen::Vector2d distorted;
    Eigen::Matrix2d J2;
    Distort(undistorted,distorted,J2);

    // Shift origin and scale
    c(0) = K_(0, 0)*distorted(0) + K_(0, 2);
    c(1) = K_(1, 1)*distorted(1) + K_(1, 2);
    Eigen::Matrix2d J3; J3.setZero();
    J3(0,0) = K_(0, 0);
    J3(1,1) = K_(1, 1);

    J = J3*J2*J1;

    return true;
  }

  bool Camera::PixelToBearing(const Eigen::Vector2d& c,Eigen::Vector3d& vec) const{
    // Shift origin and scale
    Eigen::Vector2d y;
    y(0) = (c(0) - K_(0, 2)) / K_(0, 0);
    y(1) = (c(1) - K_(1, 2)) / K_(1, 1);

    // UnDistort by optimizing
    const int max_iter = 100;
    const double tolerance = 1e-10;
    Eigen::Vector2d ybar = y; // current guess (undistorted)
    Eigen::Matrix2d J;
    Eigen::Vector2d y_tmp; // current guess (distorted)
    Eigen::Vector2d e;
    Eigen::Vector2d du;
    bool success = false;
    for (int i = 0; i < max_iter; i++) {
      Distort(ybar,y_tmp,J);
      e = y - y_tmp;
      du = (J.transpose() * J).inverse() * J.transpose() * e;
      ybar += du;
      if (e.dot(e) <= tolerance){
        success = true;
        break;
      }
    }
    if(success){
      y = ybar;
      vec = Eigen::Vector3d(y(0),y(1),1.0).normalized();
    }
    return success;
  }

  void Camera::TestCameraModel(){
//    double d = 1e-4;
//    LWF::NormalVectorElement b_s;
//    LWF::NormalVectorElement b_s1;
//    LWF::NormalVectorElement b_s2;
//    Eigen::Vector3d v_s;
//    Eigen::Vector3d v_s1;
//    Eigen::Vector3d v_s2;
//    LWF::NormalVectorElement b_e;
//    Eigen::Matrix2d J1;
//    Eigen::Matrix2d J1_FD;
//    Eigen::Matrix<double,2,3> J2;
//    Eigen::Matrix<double,2,3> J2_FD;
//    Eigen::Vector2d p_s;
//    Eigen::Vector2d p_s1;
//    Eigen::Vector2d p_s2;
//    Eigen::Vector2d p_s3;
//    Eigen::Vector2d diff;
//    for(unsigned int s = 1; s<10;){
//      b_s.setRandom(s);
//      if(b_s.getVec()(2)<0) b_s = b_s.inverted();
//      BearingToPixel(b_s,p_s,J1);
//      pixelToBearing(p_s,b_e);
//      b_s.boxMinus(b_e,diff);
//      std::cout << b_s.getVec().transpose() << std::endl;
//      std::cout << "Error after back and forward mapping: " << diff.norm() << std::endl;
//      diff = Eigen::Vector2d(d,0);
//      b_s.boxPlus(diff,b_s1);
//      BearingToPixel(b_s1,p_s1);
//      J1_FD(0,0) = static_cast<double>((p_s1-p_s).x)/d;
//      J1_FD(1,0) = static_cast<double>((p_s1-p_s).y)/d;
//      diff = Eigen::Vector2d(0,d);
//      b_s.boxPlus(diff,b_s2);
//      BearingToPixel(b_s2,p_s2);
//      J1_FD(0,1) = static_cast<double>((p_s2-p_s).x)/d;
//      J1_FD(1,1) = static_cast<double>((p_s2-p_s).y)/d;
//      std::cout << J1 << std::endl;
//      std::cout << J1_FD << std::endl;
//
//      v_s = b_s.getVec();
//      BearingToPixel(v_s,p_s,J2);
//      BearingToPixel(v_s + Eigen::Vector3d(d,0,0),p_s1);
//      BearingToPixel(v_s + Eigen::Vector3d(0,d,0),p_s2);
//      BearingToPixel(v_s + Eigen::Vector3d(0,0,d),p_s3);
//      J2_FD(0,0) = static_cast<double>((p_s1-p_s).x)/d;
//      J2_FD(1,0) = static_cast<double>((p_s1-p_s).y)/d;
//      J2_FD(0,1) = static_cast<double>((p_s2-p_s).x)/d;
//      J2_FD(1,1) = static_cast<double>((p_s2-p_s).y)/d;
//      J2_FD(0,2) = static_cast<double>((p_s3-p_s).x)/d;
//      J2_FD(1,2) = static_cast<double>((p_s3-p_s).y)/d;
//      std::cout << J2 << std::endl;
//      std::cout << J2_FD << std::endl;
//    }
  }
}
