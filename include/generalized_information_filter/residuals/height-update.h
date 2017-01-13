#ifndef GIF_HEIGHTUPDATE_HPP_
#define GIF_HEIGHTUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"
#include "generalized_information_filter/measurements/height-meas.h"

namespace GIF {

/*! \brief Pose Update with time alignment
 *         Builds a residual between current estimate pose and an external measured pose. The state
 *         is parametrized by the position (IrIB) and attitude (qIB). Additionally, the offset
 *         between the inertial frame is coestimated (IrIJ and qIJ). The extrinsics calibration of
 *         the pose measurement is also coestimated  (BrBC and qBC).
 *
 *         Coordinate frames:
 *           B: Body
 *           C: Pose measurement sensor
 *           J: Pose measurement reference inertial frame
 *           I: Estimation reference inertial frame
 */
class HeightUpdate : public UnaryUpdate<ElementPack<double>,
    ElementPack<double,Vec3>, ElementPack<double>, HeightMeas> {
 public:
  HeightUpdate(const std::string& name,
             const std::array<std::string,1>& errorName = {"z"},
             const std::array<std::string,2>& stateName = {"zRef", "IrIB"},
             const std::array<std::string,1>& noiseName = {"z"})
       : mtUnaryUpdate(name, errorName, stateName, noiseName){
  }

  virtual ~HeightUpdate() {
  }

  void Eval(double& z_inn, const double& zRef_cur, const Vec3& IrIB_cur, const double& z_noi) const {
    z_inn = IrIB_cur(2) - zRef_cur + z_noi - meas_->z_;
  }

  void JacCur(MatX& J, const double& zRef_cur, const Vec3& IrIB_cur, const double& z_noi) const {
    J.setZero();
    GetJacBlockCur<HEI, HEI>(J) = -Mat<1>::Identity();
    GetJacBlockCur<HEI, POS>(J) = Vec3(0,0,1).transpose();
  }

  void JacNoi(MatX& J, const double& zRef_cur, const Vec3& IrIB_cur, const double& z_noi) const {
    J.setZero();
    GetJacBlockNoi<HEI, HEI>(J) = Mat<1>::Identity();
  }

 protected:
  enum Elements {
    HEI,
    POS
  };
};

}

#endif /* GIF_HEIGHTUPDATE_HPP_ */
