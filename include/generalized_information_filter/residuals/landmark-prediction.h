#ifndef GIF_LANDMARKPREDICTION_HPP_
#define GIF_LANDMARKPREDICTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/prediction.h"
#include "generalized_information_filter/residuals/imu-prediction.h"

namespace GIF {

/*! \brief Imu Prediction with robocentric landmark.
 *         Derives from ImuPrediction. Add robocentric landmark which are propagated  based on
 *         the estimated rotational rate. TODO: make it dpeend on ImuPrediction
 */
template<int NumLandmark>
class RobocentricLandmarkPrediction
      : public Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
                          ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
                          ImuMeas> {
 public:
  using mtPrediction = Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
      ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
      ImuMeas>;
  using mtPrediction::meas_;
  using mtPrediction::dt_;
  RobocentricLandmarkPrediction()
    : mtPrediction({"IrIB", "BvB", "BwB_bias", "BfB_bias", "qIB", "BrBL"},
                   {"IrIB", "BvB", "BwB_bias", "BfB_bias", "qIB", "BrBL"}),
                   Ig_(0, 0, -9.81){
  }
  virtual ~RobocentricLandmarkPrediction() {
  }
  enum Elements {POS, VEL, GYB, ACB, ATT, FEA};
  void Predict(Vec3& IrIB_cur, Vec3& BvB_cur, Vec3& BwB_bias_cur, Vec3& BfB_bias_cur, Quat& qIB_cur,
               std::array<Vec3,NumLandmark>& BrBL_cur, const Vec3& IrIB_pre, const Vec3& BvB_pre,
               const Vec3& BwB_bias_pre, const Vec3& BfB_bias_pre, const Quat& qIB_pre,
               const std::array<Vec3,NumLandmark>& BrBL_pre, const Vec3& IrIB_noi,
               const Vec3& BvB_noi, const Vec3& BwB_bias_noi, const Vec3& BfB_bias_noi,
               const Vec3& qIB_noi, const std::array<Vec3,NumLandmark>& BrBL_noi) const {
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    Quat dQ = dQ.exponentialMap(dOmega);
    IrIB_cur = IrIB_pre + dt_ * (qIB_pre.rotate(BvB_pre) + IrIB_noi / sqrt(dt_));
    BvB_cur = (Mat3::Identity() - gSM(dOmega)) * BvB_pre + dt_ * (BfB + qIB_pre.inverseRotate(Ig_));
    BwB_bias_cur = BwB_bias_pre + BwB_bias_noi * sqrt(dt_);
    BfB_bias_cur = BfB_bias_pre + BfB_bias_noi * sqrt(dt_);
    qIB_cur = qIB_pre * dQ;
    for(int i=0;i<NumLandmark;i++){
      BrBL_cur[i] = (Mat3::Identity() - gSM(dOmega)) * BrBL_pre[i]
                    - dt_ * BvB_pre + BrBL_noi[i] * sqrt(dt_);
    }
  }
  void PredictJacPre(MatX& J, const Vec3& IrIB_pre, const Vec3& BvB_pre, const Vec3& BwB_bias_pre,
                              const Vec3& BfB_bias_pre, const Quat& qIB_pre, const std::array<Vec3,NumLandmark>& BrBL_pre, const Vec3& IrIB_noi,
                              const Vec3& BvB_noi, const Vec3& BwB_bias_noi,
                              const Vec3& BfB_bias_noi, const Vec3& qIB_noi, const std::array<Vec3,NumLandmark>& BrBL_noi) const {
    J.setZero();
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    this->template GetJacBlockPre<POS, POS>(J) = Mat3::Identity();
    this->template GetJacBlockPre<POS, VEL>(J) = dt_ * RotMat(qIB_pre).matrix();
    this->template GetJacBlockPre<POS, ATT>(J) = -dt_ * gSM(qIB_pre.rotate(BvB_pre));
    this->template GetJacBlockPre<VEL, VEL>(J) = (Mat3::Identity() - gSM(dOmega));
    this->template GetJacBlockPre<VEL, GYB>(J) = -dt_ * gSM(BvB_pre);
    this->template GetJacBlockPre<VEL, ACB>(J) = -dt_ * Mat3::Identity();
    this->template GetJacBlockPre<VEL, ATT>(J) = dt_ * RotMat(qIB_pre).matrix().transpose() * gSM(Ig_);
    this->template GetJacBlockPre<GYB, GYB>(J) = Mat3::Identity();
    this->template GetJacBlockPre<ACB, ACB>(J) = Mat3::Identity();
    this->template GetJacBlockPre<ATT, GYB>(J) = -dt_ * RotMat(qIB_pre).matrix() * GammaMat(dOmega);
    this->template GetJacBlockPre<ATT, ATT>(J) = Mat3::Identity();
    for(int i=0;i<NumLandmark;i++){
      this->template GetJacBlockNoi<FEA, VEL>(J).template block<3,3>(i*3,0) = -dt_ * Mat3::Identity();
      this->template GetJacBlockNoi<FEA, GYB>(J).template block<3,3>(i*3,0) = -dt_ * gSM(BrBL_pre[i]);
      this->template GetJacBlockNoi<FEA, FEA>(J).template block<3,3>(i*3,i*3) = (Mat3::Identity() - gSM(dOmega));
    }
  }
  void PredictJacNoi(MatX& J, const Vec3& IrIB_pre, const Vec3& BvB_pre, const Vec3& BwB_bias_pre,
                              const Vec3& BfB_bias_pre, const Quat& qIB_pre, const std::array<Vec3,NumLandmark>& BrBL_pre, const Vec3& IrIB_noi,
                              const Vec3& BvB_noi, const Vec3& BwB_bias_noi,
                              const Vec3& BfB_bias_noi, const Vec3& qIB_noi, const std::array<Vec3,NumLandmark>& BrBL_noi) const {
    J.setZero();
    const Vec3 BwB = meas_->BwB_ - BwB_bias_pre + qIB_noi / sqrt(dt_);
    const Vec3 BfB = meas_->BfB_ - BfB_bias_pre + BvB_noi / sqrt(dt_);
    const Vec3 dOmega = dt_ * BwB;
    this->template GetJacBlockNoi<POS, POS>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<VEL, VEL>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<VEL, ATT>(J) = sqrt(dt_) * gSM(BvB_pre);
    this->template GetJacBlockNoi<GYB, GYB>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<ACB, ACB>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<ATT, ATT>(J) = sqrt(dt_) * RotMat(qIB_pre).matrix() * GammaMat(dOmega);
    for(int i=0;i<NumLandmark;i++){
      this->template GetJacBlockNoi<FEA, ATT>(J).template block<3,3>(i*3,0) = sqrt(dt_) * gSM(BrBL_pre[i]);
      this->template GetJacBlockNoi<FEA, FEA>(J).template block<3,3>(i*3,i*3) = sqrt(dt_) * Mat3::Identity();
    }
  }

 protected:
  const Vec3 Ig_;
};

}
#endif /* GIF_LANDMARKPREDICTION_HPP_ */
