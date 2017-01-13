#ifndef GIF_LANDMARKPREDICTIONAUGMENTED_HPP_
#define GIF_LANDMARKPREDICTIONAUGMENTED_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/measurements/imu-meas.h"
#include "generalized_information_filter/prediction.h"

namespace GIF {

/*! \brief Imu Prediction with robocentric landmark.
 *         Derives from ImuPrediction. Add robocentric landmark which are propagated  based on
 *         the estimated rotational rate. TODO: make it dpeend on ImuPrediction
 */
template<int NumLandmark>
class RobocentricLandmarkPredictionAugmented
      : public BinaryResidual<
          ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
          ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
          ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
          ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
          ImuMeas> {
 public:
  typedef std::array<std::string,7> String7;
  using mtResidual = BinaryResidual<
      ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
      ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
      ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
      ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
      ImuMeas>;
  using mtResidual::meas_;
  using mtResidual::dt_;
  using mtResidual::GetJacBlockPre;
  using mtResidual::GetJacBlockCur;
  using mtResidual::GetJacBlockNoi;
  RobocentricLandmarkPredictionAugmented(const std::string& name,
            const String7& innName = {"IrIM", "MvM", "MwM", "MwM_bias", "MfM_bias", "qIM", "MrML"},
            const String7& preName = {"IrIM", "MvM", "MwM", "MwM_bias", "MfM_bias", "qIM", "MrML"},
            const String7& curName = {"IrIM", "MvM", "MwM", "MwM_bias", "MfM_bias", "qIM", "MrML"},
            const String7& noiName = {"IrIM", "MvM", "MwM", "MwM_bias", "MfM_bias", "qIM", "MrML"})
    : mtResidual(name,innName,preName,curName,noiName,false,true,true),
      Ig_(0, 0, -9.81){
    for(int i=0;i<NumLandmark;i++){
      propagation_flags_[i] = true;
    }
  }
  virtual ~RobocentricLandmarkPredictionAugmented() {
  }
  enum Elements {POS, VEL, ROR, GYB, ACB, ATT, FEA};
  void Eval(Vec3& IrIM_inn, Vec3& MvM_inn, Vec3& MwM_inn, Vec3& MwM_bias_inn, Vec3& MfM_bias_inn,
            Vec3& qIM_inn, std::array<Vec3,NumLandmark>& MrML_inn,
            const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre,
            const Vec3& MwM_bias_pre, const Vec3& MfM_bias_pre,
            const Quat& qIM_pre, const std::array<Vec3,NumLandmark>& MrML_pre,
            const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur,
            const Vec3& MwM_bias_cur, const Vec3& MfM_bias_cur,
            const Quat& qIM_cur, const std::array<Vec3,NumLandmark>& MrML_cur,
            const Vec3& IrIM_noi, const Vec3& MvM_noi, const Vec3& MwM_noi,
            const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
            const Vec3& qIM_noi, const std::array<Vec3,NumLandmark>& MrML_noi) const {
    ComputeBiasCorrectedAccMeas(MfM_bias_pre, MvM_noi);
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur + qIM_noi * sqrt(dt_));
    IrIM_inn = IrIM_pre + dt_ * (qIM_pre.rotate(MvM_pre) + IrIM_noi / sqrt(dt_)) - IrIM_cur;
    MvM_inn = (Mat3::Identity() - gSM(dt_ * MwM_cur)) * MvM_pre
              + dt_ * (MfM_cor_ + qIM_pre.inverseRotate(Ig_)) - MvM_cur;
    MwM_bias_inn = MwM_bias_pre + MwM_bias_noi * sqrt(dt_) - MwM_bias_cur;
    MfM_bias_inn = MfM_bias_pre + MfM_bias_noi * sqrt(dt_) - MfM_bias_cur;
    qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        MrML_inn[i] = (Mat3::Identity() - gSM(dt_ * MwM_cur)) * MrML_pre[i]
                      - dt_ * MvM_pre + MrML_noi[i] * sqrt(dt_) - MrML_cur[i];
      } else {
        MrML_inn[i] = MrML_noi[i] * sqrt(dt_) -  MrML_cur[i];
      }
    }
    MwM_inn = meas_->MwM_ - (MwM_cur + MwM_bias_cur + MwM_noi / sqrt(dt_));
  }
  void JacPre(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre,
              const Vec3& MwM_bias_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur,
              const Vec3& MwM_bias_cur, const Vec3& MfM_bias_cur,
              const Quat& qIM_cur, const std::array<Vec3,NumLandmark>& MrML_cur,
              const Vec3& IrIM_noi, const Vec3& MvM_noi, const Vec3& MwM_noi,
              const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
              const Vec3& qIM_noi, const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    ComputeBiasCorrectedAccMeas(MfM_bias_pre, MvM_noi);
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockPre<POS, POS>(J) = Mat3::Identity();
    this->template GetJacBlockPre<POS, VEL>(J) = dt_ * RotMat(qIM_pre).matrix();
    this->template GetJacBlockPre<POS, ATT>(J) = -dt_ * gSM(qIM_pre.rotate(MvM_pre));
    this->template GetJacBlockPre<VEL, VEL>(J) = (Mat3::Identity() - gSM(dt_ * MwM_cur));
    this->template GetJacBlockPre<VEL, ACB>(J) = -dt_ * Mat3::Identity();
    this->template GetJacBlockPre<VEL, ATT>(J) = dt_ * RotMat(qIM_pre).matrix().transpose() *
        gSM(Ig_);
    this->template GetJacBlockPre<GYB, GYB>(J) = Mat3::Identity();
    this->template GetJacBlockPre<ACB, ACB>(J) = Mat3::Identity();
    this->template GetJacBlockPre<ATT, ATT>(J) = GammaMat(qIM_inn).inverse();
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockPre<FEA, VEL>(J).template block<3,3>(i*3,0) =
            -dt_ * Mat3::Identity();
        this->template GetJacBlockPre<FEA, FEA>(J).template block<3,3>(i*3,i*3) =
            (Mat3::Identity() - gSM(dt_ * MwM_cur));
      }
    }
  }
  void JacCur(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre,
              const Vec3& MwM_bias_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur,
              const Vec3& MwM_bias_cur, const Vec3& MfM_bias_cur,
              const Quat& qIM_cur, const std::array<Vec3,NumLandmark>& MrML_cur,
              const Vec3& IrIM_noi, const Vec3& MvM_noi, const Vec3& MwM_noi,
              const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
              const Vec3& qIM_noi, const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    ComputeBiasCorrectedAccMeas(MfM_bias_pre, MvM_noi);
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockCur<POS, POS>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<VEL, VEL>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<VEL, ROR>(J) = dt_ * gSM(MvM_pre);
    this->template GetJacBlockCur<GYB, GYB>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<ACB, ACB>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<ATT, ATT>(J) = -GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre * dQ * qIM_cur.inverted()).matrix();
    this->template GetJacBlockCur<ATT, ROR>(J) = dt_ * GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cur);
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockCur<FEA, FEA>(J).template block<3,3>(i*3,i*3) = -Mat3::Identity();
        this->template GetJacBlockCur<FEA, ROR>(J).template block<3,3>(i*3,0) =
            dt_ * gSM(MrML_pre[i]);
      }
    }
    this->template GetJacBlockCur<ROR, ROR>(J) = -Mat3::Identity();
    this->template GetJacBlockCur<ROR, GYB>(J) = -Mat3::Identity();
  }
  void JacNoi(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre,
              const Vec3& MwM_bias_pre, const Vec3& MfM_bias_pre,
              const Quat& qIM_pre, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur,
              const Vec3& MwM_bias_cur, const Vec3& MfM_bias_cur,
              const Quat& qIM_cur, const std::array<Vec3,NumLandmark>& MrML_cur,
              const Vec3& IrIM_noi, const Vec3& MvM_noi, const Vec3& MwM_noi,
              const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
              const Vec3& qIM_noi, const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    ComputeBiasCorrectedAccMeas(MfM_bias_pre, MvM_noi);
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockNoi<POS, POS>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<VEL, VEL>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<GYB, GYB>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<ACB, ACB>(J) = sqrt(dt_) * Mat3::Identity();
    this->template GetJacBlockNoi<ATT, ATT>(J) = sqrt(dt_) * GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cur);
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
      }
      this->template GetJacBlockNoi<FEA, FEA>(J).template block<3,3>(i*3,i*3) =
          sqrt(dt_) * Mat3::Identity();
    }
    this->template GetJacBlockNoi<ROR, ROR>(J) = -1/sqrt(dt_) * Mat3::Identity();
  }
  void ComputeBiasCorrectedAccMeas(const Vec3& MfM_bias_pre, const Vec3& MvM_noi) const{
    MfM_cor_ = meas_->MfM_ - MfM_bias_pre + MvM_noi / sqrt(dt_);
  }
  Vec3 GetMeasuredMwM() const{
    return meas_->MwM_;
  }
  void SetPropagationFlag(int i, bool flag) {
    propagation_flags_[i] = flag;
  }

 protected:
  std::array<bool,NumLandmark> propagation_flags_;
  const Vec3 Ig_;
  mutable Vec3 MfM_cor_;
};

}
#endif /* GIF_LANDMARKPREDICTIONAUGMENTED_HPP_ */
