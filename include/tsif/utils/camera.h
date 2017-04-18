#ifndef TSIF_CAMERA_H_
#define TSIF_CAMERA_H_

#include "tsif/utils/common.h"

namespace tsif{

class Camera{
 public:
  enum ModelType{
    RADTAN,
    EQUIDIST
  } type_;

  Eigen::Matrix3d K_;

  double k1_,k2_,k3_,k4_,k5_,k6_;
  double p1_,p2_,s1_,s2_,s3_,s4_;

  Camera();
  ~Camera();
  void LoadCameraMatrix(const std::string& filename);
  void LoadRadtan(const std::string& filename);
  void LoadEquidist(const std::string& filename);
  void Load(const std::string& filename);
  void DistortRadtan(const Eigen::Vector2d& in, Eigen::Vector2d& out) const;
  void DistortRadtan(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const;
  void DistortEquidist(const Eigen::Vector2d& in, Eigen::Vector2d& out) const;
  void DistortEquidist(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const;
  void Distort(const Eigen::Vector2d& in, Eigen::Vector2d& out) const;
  void Distort(const Eigen::Vector2d& in, Eigen::Vector2d& out, Eigen::Matrix2d& J) const;
  bool BearingToPixel(const Eigen::Vector3d& vec, Eigen::Vector2d& c) const;
  bool BearingToPixel(const Eigen::Vector3d& vec, Eigen::Vector2d& c, Eigen::Matrix<double,2,3>& J) const;
  bool PixelToBearing(const Eigen::Vector2d& c,Eigen::Vector3d& vec) const;
  void TestCameraModel();
};

}


#endif /* TSIF_CAMERA_H_ */
