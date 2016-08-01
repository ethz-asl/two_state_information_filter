/*
 * Common.hpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Bloeschm
 */

#ifndef GIF_COMMON_HPP_
#define GIF_COMMON_HPP_

#include <Eigen/Dense>
#include "kindr/Core"
#include <iostream>

typedef kindr::RotationQuaternionPD QPD;
typedef kindr::RotationMatrixPD MPD;
typedef Eigen::Vector3d V3D;
typedef Eigen::Matrix3d M3D;
typedef Eigen::VectorXd VXD;
typedef Eigen::MatrixXd MXD;
inline M3D gSM(const V3D& vec){
  return kindr::getSkewMatrixFromVector(vec);
}

static void enforceSymmetry(MXD& mat){
  mat = 0.5*(mat+mat.transpose()).eval();
}

inline M3D Lmat (const V3D& a) {
  return kindr::getJacobianOfExponentialMap(a);
}

namespace GIF{
  enum FilteringMode{
    ModeEKF,
    ModeUKF,
    ModeIEKF
  };
}

#endif /* GIF_COMMON_HPP_ */
