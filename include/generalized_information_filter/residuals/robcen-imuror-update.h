#ifndef GIF_ROBCENIMURORUPDATE_HPP_
#define GIF_ROBCENIMURORUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/residuals/ror-meas.h"
#include "generalized_information_filter/unary-update.h"

namespace GIF {

/*! \brief Imu based robocentric rotational rate update
 */
class RobcenImurorUpdate : public BinaryResidual<ElementPack<Vec3>, ElementPack<>,
                                        ElementPack<Vec3,Vec3>, ElementPack<Vec3>, RorMeas> {
 public:
  using mtResidual = BinaryResidual<ElementPack<Vec3>, ElementPack<>,
      ElementPack<Vec3,Vec3>, ElementPack<Vec3>, RorMeas>;
  RobcenImurorUpdate(const std::string& name,
             const std::array<std::string,1>& innName = {"MwM"},
             const std::array<std::string,0>& preName = {},
             const std::array<std::string,2>& curName = {"MwM", "MwM_bias"},
             const std::array<std::string,1>& noiName = {"MwM"})
       : mtResidual(name, innName, preName, curName, noiName,false,true,true){
  }
  virtual ~RobcenImurorUpdate() {
  }
  enum Elements {ROR, GYB};
  void Eval(Vec3& MwM_inn, const Vec3& MwM_cur, const Vec3& MwM_bias_cur,
            const Vec3& MwM_noi) const {
    MwM_inn = meas_->MwM_ - (MwM_cur + MwM_bias_cur + MwM_noi / sqrt(dt_));
  }

  void JacPre(MatX& J, const Vec3& MwM_cur, const Vec3& MwM_bias_cur,
              const Vec3& MwM_noi) const {
  }

  void JacCur(MatX& J, const Vec3& MwM_cur, const Vec3& MwM_bias_cur,
              const Vec3& MwM_noi) const {
    J.setZero();
    this->template GetJacBlockCur<ROR, ROR>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<ROR, GYB>(J) = -Mat3::Identity();
  }

  void JacNoi(MatX& J, const Vec3& MwM_cur, const Vec3& MwM_bias_cur,
              const Vec3& MwM_noi) const {
    J.setZero();
    this->template GetJacBlockNoi<ROR, ROR>(J) = -1/sqrt(dt_) * Mat3::Identity();
  }
};

}

#endif /* GIF_ROBCENIMURORUPDATE_HPP_ */
