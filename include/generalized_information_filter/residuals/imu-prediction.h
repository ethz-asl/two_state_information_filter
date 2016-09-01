#ifndef GIF_IMUPREDICTION_HPP_
#define GIF_IMUPREDICTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/prediction.h"

namespace GIF {

/*! \brief IMU Measurement
 *         ElementVector that can be used to hold IMU measurements.
 */
class IMUMeas : public ElementVector {
 public:
  IMUMeas(const Vec3& BwB = Vec3(0, 0, 0), const Vec3& BfB = Vec3(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Vec3>({"BwB", "BfB"}))),
        BwB_(ElementVector::GetValue<Vec3>("BwB")),
        BfB_(ElementVector::GetValue<Vec3>("BfB")) {
    BwB_ = BwB;
    BfB_ = BfB;
  }
  Vec3& BwB_;
  Vec3& BfB_;
};

/*! \brief IMU Prediction.
 *         Predicts the current state based on the previous state and the IMU measurement.
 *         Coordinate frames: I (world), B (IMU).
 *         Chosen state parametrisation: Position (IrIB), Velocity (BvB), Attitude (qIB),
 *                                       Gyroscope bias (BwB_bias), Accelerometer bias (BfB_bias).
 */
class IMUPrediction : public Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat>,
                                        ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3>,
                                        IMUMeas> {
 public:
  IMUPrediction()
      : mtPrediction({"IrIB", "BvB", "BwB_bias", "BfB_bias", "qIB"},
                     {"IrIB", "BvB", "BwB_bias", "BfB_bias", "qIB"}),
        Ig_(0, 0, -9.81),
        dt_(0.1){
  }
  virtual ~IMUPrediction() {
  }
  void Predict(Vec3& IrIB_cur, Vec3& BvB_cur, Vec3& BwB_bias_cur, Vec3& BfB_bias_cur, Quat& qIB_cur,
               const Vec3& IrIB_pre, const Vec3& BvB_pre, const Vec3& BwB_bias_pre,
               const Vec3& BfB_bias_pre, const Quat& qIB_pre, const Vec3& IrIB_noi,
               const Vec3& BvB_noi, const Vec3& BwB_bias_noi, const Vec3& BfB_bias_noi,
               const Vec3& qIB_noi) const {
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    Quat dQ = dQ.exponentialMap(dOmega);
    IrIB_cur = IrIB_pre + dt_ * (qIB_pre.rotate(BvB_pre) + IrIB_noi / sqrt(dt_));
    BvB_cur = (Mat3::Identity() - gSM(dOmega)) * BvB_pre + dt_ * (BfB + qIB_pre.inverseRotate(Ig_));
    BwB_bias_cur = BwB_bias_pre + BwB_bias_noi * sqrt(dt_);
    BfB_bias_cur = BfB_bias_pre + BfB_bias_noi * sqrt(dt_);
    qIB_cur = qIB_pre * dQ;
  }
  void PredictJacPre(MatX& J, const Vec3& IrIB_pre, const Vec3& BvB_pre, const Vec3& BwB_bias_pre,
                              const Vec3& BfB_bias_pre, const Quat& qIB_pre, const Vec3& IrIB_noi,
                              const Vec3& BvB_noi, const Vec3& BwB_bias_noi,
                              const Vec3& BfB_bias_noi, const Vec3& qIB_noi) const {
    J.setZero();
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    SetJacBlockPre<POS, POS>(J, Mat3::Identity());
    SetJacBlockPre<POS, VEL>(J, dt_ * RotMat(qIB_pre).matrix());
    SetJacBlockPre<POS, ATT>(J, -dt_ * gSM(qIB_pre.rotate(BvB_pre)));
    SetJacBlockPre<VEL, VEL>(J, (Mat3::Identity() - gSM(dOmega)));
    SetJacBlockPre<VEL, GYB>(J, -dt_ * gSM(BvB_pre));
    SetJacBlockPre<VEL, ACB>(J, -dt_ * Mat3::Identity());
    SetJacBlockPre<VEL, ATT>(J, dt_ * RotMat(qIB_pre).matrix().transpose() * gSM(Ig_));
    SetJacBlockPre<GYB, GYB>(J, Mat3::Identity());
    SetJacBlockPre<ACB, ACB>(J, Mat3::Identity());
    SetJacBlockPre<ATT, GYB>(J, -dt_ * RotMat(qIB_pre).matrix() * GammaMat(dOmega));
    SetJacBlockPre<ATT, ATT>(J, Mat3::Identity());
  }
  void PredictJacNoi(MatX& J, const Vec3& IrIB_pre, const Vec3& BvB_pre, const Vec3& BwB_bias_pre,
                              const Vec3& BfB_bias_pre, const Quat& qIB_pre, const Vec3& IrIB_noi,
                              const Vec3& BvB_noi, const Vec3& BwB_bias_noi,
                              const Vec3& BfB_bias_noi, const Vec3& qIB_noi) const {
    J.setZero();
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    SetJacBlockNoi<POS, POS>(J, sqrt(dt_) * Mat3::Identity());
    SetJacBlockNoi<VEL, VEL>(J, sqrt(dt_) * Mat3::Identity());
    SetJacBlockNoi<VEL, ATT>(J, sqrt(dt_) * gSM(BvB_pre));
    SetJacBlockNoi<GYB, GYB>(J, sqrt(dt_) * Mat3::Identity());
    SetJacBlockNoi<ACB, ACB>(J, sqrt(dt_) * Mat3::Identity());
    SetJacBlockNoi<ATT, ATT>(J, sqrt(dt_) * RotMat(qIB_pre).matrix() * GammaMat(dOmega));
  }

 protected:
  double dt_;
  const Vec3 Ig_;
  enum Elements {POS, VEL, GYB, ACB, ATT};
};

}
#endif /* GIF_IMUPREDICTION_HPP_ */
