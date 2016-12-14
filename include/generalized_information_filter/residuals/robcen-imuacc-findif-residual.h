#ifndef GIF_ROBCENIMUACCFINDIFRESIDUAL_HPP_
#define GIF_ROBCENIMUACCFINDIFRESIDUAL_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/residuals/acc-meas.h"
#include "generalized_information_filter/binary-residual.h"

namespace GIF {

/*! \brief Imu based robocentric acceleration finite difference residual
 */
class RobcenImuaccFindifResidual
      : public BinaryResidual<ElementPack<Vec3>, ElementPack<Vec3, Vec3, Quat>,
                              ElementPack<Vec3, Vec3>, ElementPack<Vec3>, AccMeas> {
 public:
  using mtResidual = BinaryResidual<ElementPack<Vec3>, ElementPack<Vec3, Vec3, Quat>,
                                    ElementPack<Vec3, Vec3>, ElementPack<Vec3>, AccMeas>;
  RobcenImuaccFindifResidual(const std::string& name,
            const std::array<std::string,1>& innName = {"MvM"},
            const std::array<std::string,3>& preName = {"MvM", "MfM_bias", "qIM"},
            const std::array<std::string,2>& curName = {"MvM", "MwM"},
            const std::array<std::string,1>& noiName = {"MvM"})
    : mtResidual(name,innName,preName,curName,noiName,false,true,true),
      Ig_(0, 0, -9.81){
  }
  virtual ~RobcenImuaccFindifResidual() {
  }
  enum Elements {VEL, ROR, ACB, ATT};
  enum ElementsInn {VEL_INN};
  enum ElementsPre {VEL_PRE, ACB_PRE, ATT_PRE};
  enum ElementsCur {VEL_CUR, ROR_CUR};
  enum ElementsNoi {VEL_NOI};
  void Eval(Vec3& MvM_inn, const Vec3& MvM_pre, const Vec3& MfM_bias_pre,
            const Quat& qIM_pre, const Vec3& MvM_cur, const Vec3& MwM_cur, const Vec3& MvM_noi) const {
    Vec3 MfM_cor = meas_->MfM_ - MfM_bias_pre + MvM_noi / sqrt(dt_);
    MvM_inn = (Mat3::Identity() - gSM(dt_ * MwM_cur)) * MvM_pre
              + dt_ * (MfM_cor + qIM_pre.inverseRotate(Ig_)) - MvM_cur;
  }
  void JacPre(MatX& J, const Vec3& MvM_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const Vec3& MvM_cur, const Vec3& MwM_cur, const Vec3& MvM_noi) const {
    J.setZero();
    this->template GetJacBlockPre<VEL_INN, VEL_PRE>(J) = (Mat3::Identity() - gSM(dt_ * MwM_cur));
    this->template GetJacBlockPre<VEL_INN, ACB_PRE>(J) = -dt_ * Mat3::Identity();
    this->template GetJacBlockPre<VEL_INN, ATT_PRE>(J) = dt_ * RotMat(qIM_pre).matrix().transpose() *
        gSM(Ig_);
  }
  void JacCur(MatX& J, const Vec3& MvM_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const Vec3& MvM_cur, const Vec3& MwM_cur, const Vec3& MvM_noi) const {
    J.setZero();
    this->template GetJacBlockCur<VEL_INN, VEL_CUR>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<VEL_INN, ROR_CUR>(J) = dt_ * gSM(MvM_pre);
  }
  void JacNoi(MatX& J, const Vec3& MvM_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const Vec3& MvM_cur, const Vec3& MwM_cur, const Vec3& MvM_noi) const {
    J.setZero();
    this->template GetJacBlockNoi<VEL_INN, VEL_NOI>(J) = sqrt(dt_) * Mat3::Identity();
  }

 protected:
  const Vec3 Ig_;
};

}
#endif /* GIF_ROBCENIMUACCFINDIFRESIDUAL_HPP_ */
