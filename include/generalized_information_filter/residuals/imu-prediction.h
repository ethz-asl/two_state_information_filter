#ifndef GIF_IMUPREDICTION_HPP_
#define GIF_IMUPREDICTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/prediction.h"

namespace GIF {

class IMUMeas : public ElementVector {
 public:
  IMUMeas(const Vec3& gyr = Vec3(0, 0, 0), const Vec3& acc = Vec3(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Vec3>({"gyr", "acc"}))),
        gyr_(ElementVector::GetValue<Vec3>("gyr")),
        acc_(ElementVector::GetValue<Vec3>("acc")) {
    gyr_ = gyr;
    acc_ = acc;
  }
  Vec3& gyr_;
  Vec3& acc_;
};

class IMUPrediction : public Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat>,
                                        ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3>,
                                        IMUMeas> {
 public:
  IMUPrediction()
      : mtPrediction({"pos", "vel", "gyb", "acb", "att"},
                     {"pos", "vel", "gyb", "acb", "att"}),
        g_(0, 0, -9.81),
        dt_(0.1){
  }
  virtual ~IMUPrediction() {
  }
  void predict(Vec3& posCur, Vec3& velCur, Vec3& gybCur, Vec3& acbCur, Quat& attCur,
               const Vec3& posPre, const Vec3& velPre, const Vec3& gybPre,
               const Vec3& acbPre, const Quat& attPre, const Vec3& posNoi,
               const Vec3& velNoi, const Vec3& gybNoi, const Vec3& acbNoi,
               const Vec3& attNoi) const {
    const Vec3 gyr = meas_->gyr_ - gybPre + attNoi / sqrt(dt_);
    const Vec3 acc = meas_->acc_ - acbPre + velNoi / sqrt(dt_);
    const Vec3 dOmega = dt_ * gyr;
    Quat dQ = dQ.exponentialMap(dOmega);
    posCur = posPre + dt_ * (attPre.rotate(velPre) + posNoi / sqrt(dt_));
    velCur = (Mat3::Identity() - gSM(dOmega)) * velPre + dt_ * (acc + attPre.inverseRotate(g_));
    gybCur = gybPre + gybNoi * sqrt(dt_);
    acbCur = acbPre + acbNoi * sqrt(dt_);
    attCur = attPre * dQ;
  }
  void predictJacPre(MatX& J, const Vec3& posPre, const Vec3& velPre, const Vec3& gybPre,
                              const Vec3& acbPre, const Quat& attPre, const Vec3& posNoi,
                              const Vec3& velNoi, const Vec3& gybNoi, const Vec3& acbNoi,
                              const Vec3& attNoi) const {
    J.setZero();
    const Vec3 gyr = meas_->gyr_ - gybPre + attNoi / sqrt(dt_);
    const Vec3 acc = meas_->acc_ - acbPre + velNoi / sqrt(dt_);
    const Vec3 dOmega = dt_ * gyr;
    setJacBlockPre<POS, POS>(J, Mat3::Identity());
    setJacBlockPre<POS, VEL>(J, dt_ * RotMat(attPre).matrix());
    setJacBlockPre<POS, ATT>(J, -dt_ * gSM(attPre.rotate(velPre)));
    setJacBlockPre<VEL, VEL>(J, (Mat3::Identity() - gSM(dOmega)));
    setJacBlockPre<VEL, GYB>(J, -dt_ * gSM(velPre));
    setJacBlockPre<VEL, ACB>(J, -dt_ * Mat3::Identity());
    setJacBlockPre<VEL, ATT>(J, dt_ * RotMat(attPre).matrix().transpose() * gSM(g_));
    setJacBlockPre<GYB, GYB>(J, Mat3::Identity());
    setJacBlockPre<ACB, ACB>(J, Mat3::Identity());
    setJacBlockPre<ATT, GYB>(J, -dt_ * RotMat(attPre).matrix() * GammaMat(dOmega));
    setJacBlockPre<ATT, ATT>(J, Mat3::Identity());
  }
  void predictJacNoi(MatX& J, const Vec3& posPre, const Vec3& velPre, const Vec3& gybPre,
                              const Vec3& acbPre, const Quat& attPre, const Vec3& posNoi,
                              const Vec3& velNoi, const Vec3& gybNoi, const Vec3& acbNoi,
                              const Vec3& attNoi) const {
    J.setZero();
    const Vec3 gyr = meas_->gyr_ - gybPre + attNoi / sqrt(dt_);
    const Vec3 acc = meas_->acc_ - acbPre + velNoi / sqrt(dt_);
    const Vec3 dOmega = dt_ * gyr;
    setJacBlockNoi<POS, POS>(J, sqrt(dt_) * Mat3::Identity());
    setJacBlockNoi<VEL, VEL>(J, sqrt(dt_) * Mat3::Identity());
    setJacBlockNoi<VEL, ATT>(J, sqrt(dt_) * gSM(velPre));
    setJacBlockNoi<GYB, GYB>(J, sqrt(dt_) * Mat3::Identity());
    setJacBlockNoi<ACB, ACB>(J, sqrt(dt_) * Mat3::Identity());
    setJacBlockNoi<ATT, ATT>(J, sqrt(dt_) * RotMat(attPre).matrix() * GammaMat(dOmega));
  }

 protected:
  double dt_;
  const Vec3 g_;
  enum Elements {POS, VEL, GYB, ACB, ATT};
};

}
#endif /* GIF_IMUPREDICTION_HPP_ */
