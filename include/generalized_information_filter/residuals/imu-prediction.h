#ifndef GIF_IMUPREDICTION_HPP_
#define GIF_IMUPREDICTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/measurements/imu-meas.h"
#include "generalized_information_filter/prediction.h"

namespace GIF {

/*! \brief Imu Prediction.
 *         Predicts the current state based on the previous state and the Imu measurement.
 *         Coordinate frames: I (world), M (Imu).
 *         Chosen state parametrisation: Position (IrIM), Velocity (MvM), Attitude (qIM),
 *                                       Gyroscope bias (MwM_bias), Accelerometer bias (MfM_bias).
 */
class ImuPrediction : public Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat>,
                                        ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3>,
                                        ImuMeas> {
 public:
  typedef std::array<std::string,5> String5;
  ImuPrediction(const std::string& name,
                const String5& staName = {"IrIM", "MvM", "MwM_bias", "MfM_bias", "qIM"},
                const String5& noiName = {"IrIM", "MvM", "MwM_bias", "MfM_bias", "qIM"})
      : mtPrediction(name, staName, noiName),
        Ig_(0, 0, -9.81){
  }
  virtual ~ImuPrediction() {
  }
  enum Elements {POS, VEL, GYB, ACB, ATT};
  void Predict(Vec3& IrIM_cur, Vec3& MvM_cur, Vec3& MwM_bias_cur, Vec3& MfM_bias_cur, Quat& qIM_cur,
               const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_bias_pre,
               const Vec3& MfM_bias_pre, const Quat& qIM_pre, const Vec3& IrIM_noi,
               const Vec3& MvM_noi, const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
               const Vec3& qIM_noi) const {
    ComputeBiasCorrectedImuMeas(MwM_bias_pre, qIM_noi, MfM_bias_pre, MvM_noi);
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cor_);
    IrIM_cur = IrIM_pre + dt_ * (qIM_pre.rotate(MvM_pre) + IrIM_noi / sqrt(dt_));
    MvM_cur = (Mat3::Identity() - gSM(dt_ * MwM_cor_)) * MvM_pre
              + dt_ * (MfM_cor_ + qIM_pre.inverseRotate(Ig_));
    MwM_bias_cur = MwM_bias_pre + MwM_bias_noi * sqrt(dt_);
    MfM_bias_cur = MfM_bias_pre + MfM_bias_noi * sqrt(dt_);
    qIM_cur = qIM_pre * dQ;
  }
  void PredictJacPre(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_bias_pre,
                              const Vec3& MfM_bias_pre, const Quat& qIM_pre, const Vec3& IrIM_noi,
                              const Vec3& MvM_noi, const Vec3& MwM_bias_noi,
                              const Vec3& MfM_bias_noi, const Vec3& qIM_noi) const {
    J.setZero();
    ComputeBiasCorrectedImuMeas(MwM_bias_pre, qIM_noi, MfM_bias_pre, MvM_noi);
    GetJacBlockPre<POS, POS>(J) = Mat3::Identity();
    GetJacBlockPre<POS, VEL>(J) = dt_ * RotMat(qIM_pre).matrix();
    GetJacBlockPre<POS, ATT>(J) = -dt_ * gSM(qIM_pre.rotate(MvM_pre));
    GetJacBlockPre<VEL, VEL>(J) = (Mat3::Identity() - gSM(dt_ * MwM_cor_));
    GetJacBlockPre<VEL, GYB>(J) = -dt_ * gSM(MvM_pre);
    GetJacBlockPre<VEL, ACB>(J) = -dt_ * Mat3::Identity();
    GetJacBlockPre<VEL, ATT>(J) = dt_ * RotMat(qIM_pre).matrix().transpose() * gSM(Ig_);
    GetJacBlockPre<GYB, GYB>(J) = Mat3::Identity();
    GetJacBlockPre<ACB, ACB>(J) = Mat3::Identity();
    GetJacBlockPre<ATT, GYB>(J) = -dt_ * RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cor_);
    GetJacBlockPre<ATT, ATT>(J) = Mat3::Identity();
  }
  void PredictJacNoi(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_bias_pre,
                              const Vec3& MfM_bias_pre, const Quat& qIM_pre, const Vec3& IrIM_noi,
                              const Vec3& MvM_noi, const Vec3& MwM_bias_noi,
                              const Vec3& MfM_bias_noi, const Vec3& qIM_noi) const {
    J.setZero();
    ComputeBiasCorrectedImuMeas(MwM_bias_pre, qIM_noi, MfM_bias_pre, MvM_noi);
    GetJacBlockNoi<POS, POS>(J) = sqrt(dt_) * Mat3::Identity();
    GetJacBlockNoi<VEL, VEL>(J) = sqrt(dt_) * Mat3::Identity();
    GetJacBlockNoi<VEL, ATT>(J) = sqrt(dt_) * gSM(MvM_pre);
    GetJacBlockNoi<GYB, GYB>(J) = sqrt(dt_) * Mat3::Identity();
    GetJacBlockNoi<ACB, ACB>(J) = sqrt(dt_) * Mat3::Identity();
    GetJacBlockNoi<ATT, ATT>(J) = sqrt(dt_) * RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cor_);
  }
  void ComputeBiasCorrectedImuMeas(const Vec3& MwM_bias_pre, const Vec3& qIM_noi,
                                    const Vec3& MfM_bias_pre, const Vec3& MvM_noi) const{
    MwM_cor_ = meas_->MwM_ - MwM_bias_pre + qIM_noi / sqrt(dt_);
    MfM_cor_ = meas_->MfM_ - MfM_bias_pre + MvM_noi / sqrt(dt_);
  }
  Vec3 GetMwMCor() const{
    return MwM_cor_;
  }

 protected:
  const Vec3 Ig_;
  mutable Vec3 MwM_cor_;
  mutable Vec3 MfM_cor_;
};

}
#endif /* GIF_IMUPREDICTION_HPP_ */
