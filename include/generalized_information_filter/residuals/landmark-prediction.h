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
  typedef std::array<std::string,6> String6;
  using mtPrediction =
      Prediction<ElementPack<Vec3, Vec3, Vec3, Vec3, Quat, std::array<Vec3,NumLandmark>>,
                 ElementPack<Vec3, Vec3, Vec3, Vec3, Vec3, std::array<Vec3,NumLandmark>>,
                 ImuMeas>;
  using mtPrediction::meas_;
  using mtPrediction::dt_;
  RobocentricLandmarkPrediction(
      const String6& staName = {"IrIM", "MvM", "MwM_bias", "MfM_bias", "qIM", "MrML"},
      const String6& noiName = {"IrIM", "MvM", "MwM_bias", "MfM_bias", "qIM", "MrML"})
    : mtPrediction(staName, noiName),
      imuPrediction_({staName[0], staName[1], staName[2], staName[3], staName[4]},
                     {noiName[0], noiName[1], noiName[2], noiName[3], noiName[4]}){
    for(int i=0;i<NumLandmark;i++){
      propagation_flags_[i] = true;
    }
  }
  virtual ~RobocentricLandmarkPrediction() {
  }
  enum Elements {POS, VEL, GYB, ACB, ATT, FEA};
  void Predict(Vec3& IrIM_cur, Vec3& MvM_cur, Vec3& MwM_bias_cur, Vec3& MfM_bias_cur, Quat& qIM_cur,
               std::array<Vec3,NumLandmark>& MrML_cur, const Vec3& IrIM_pre, const Vec3& MvM_pre,
               const Vec3& MwM_bias_pre, const Vec3& MfM_bias_pre, const Quat& qIM_pre,
               const std::array<Vec3,NumLandmark>& MrML_pre, const Vec3& IrIM_noi,
               const Vec3& MvM_noi, const Vec3& MwM_bias_noi, const Vec3& MfM_bias_noi,
               const Vec3& qIM_noi, const std::array<Vec3,NumLandmark>& MrML_noi) const {
    imuPrediction_.SetDt(dt_);
    imuPrediction_.SetMeas(meas_);
    imuPrediction_.Predict(IrIM_cur, MvM_cur, MwM_bias_cur, MfM_bias_cur, qIM_cur,
                           IrIM_pre, MvM_pre, MwM_bias_pre, MfM_bias_pre, qIM_pre,
                           IrIM_noi, MvM_noi, MwM_bias_noi, MfM_bias_noi, qIM_noi);
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        MrML_cur[i] = (Mat3::Identity() - gSM(dt_ * imuPrediction_.GetMwMCor())) * MrML_pre[i]
                      - dt_ * MvM_pre + MrML_noi[i] * sqrt(dt_);
      } else {
        MrML_cur[i] = MrML_noi[i] * sqrt(dt_);
      }
    }
  }
  void PredictJacPre(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_bias_pre,
                              const Vec3& MfM_bias_pre, const Quat& qIM_pre,
                              const std::array<Vec3,NumLandmark>& MrML_pre, const Vec3& IrIM_noi,
                              const Vec3& MvM_noi, const Vec3& MwM_bias_noi,
                              const Vec3& MfM_bias_noi, const Vec3& qIM_noi,
                              const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    MatX JImu(imuPrediction_.PreDefinition()->GetDim(), imuPrediction_.CurDefinition()->GetDim());
    imuPrediction_.SetDt(dt_);
    imuPrediction_.SetMeas(meas_);
    imuPrediction_.PredictJacPre(JImu,
                                 IrIM_pre, MvM_pre, MwM_bias_pre, MfM_bias_pre, qIM_pre,
                                 IrIM_noi, MvM_noi, MwM_bias_noi, MfM_bias_noi, qIM_noi);
    J.block(0,0,imuPrediction_.PreDefinition()->GetDim(), imuPrediction_.CurDefinition()->GetDim())
        = JImu;
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockNoi<FEA, VEL>(J).template block<3,3>(i*3,0) =
            -dt_ * Mat3::Identity();
        this->template GetJacBlockNoi<FEA, GYB>(J).template block<3,3>(i*3,0) =
            -dt_ * gSM(MrML_pre[i]);
        this->template GetJacBlockNoi<FEA, FEA>(J).template block<3,3>(i*3,i*3) =
            (Mat3::Identity() - gSM(dt_ * imuPrediction_.GetMwMCor()));
      }
    }
  }
  void PredictJacNoi(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_bias_pre,
                              const Vec3& MfM_bias_pre, const Quat& qIM_pre,
                              const std::array<Vec3,NumLandmark>& MrML_pre, const Vec3& IrIM_noi,
                              const Vec3& MvM_noi, const Vec3& MwM_bias_noi,
                              const Vec3& MfM_bias_noi, const Vec3& qIM_noi,
                              const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    MatX JImu(imuPrediction_.PreDefinition()->GetDim(), imuPrediction_.NoiDefinition()->GetDim());
    imuPrediction_.SetDt(dt_);
    imuPrediction_.SetMeas(meas_);
    imuPrediction_.PredictJacNoi(JImu,
                                 IrIM_pre, MvM_pre, MwM_bias_pre, MfM_bias_pre, qIM_pre,
                                 IrIM_noi, MvM_noi, MwM_bias_noi, MfM_bias_noi, qIM_noi);
    J.block(0,0,imuPrediction_.PreDefinition()->GetDim(), imuPrediction_.NoiDefinition()->GetDim())
        = JImu;
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockNoi<FEA, ATT>(J).template block<3,3>(i*3,0) =
            sqrt(dt_) * gSM(MrML_pre[i]);
      }
      this->template GetJacBlockNoi<FEA, FEA>(J).template block<3,3>(i*3,i*3) =
          sqrt(dt_) * Mat3::Identity();
    }
  }
  Vec3 GetMeasuredMwM() const{
    return meas_->MwM_;
  }

  void SetPropagationFlag(int i, bool flag) {
    propagation_flags_[i] = flag;
  }

 protected:
  std::array<bool,NumLandmark> propagation_flags_;
  mutable ImuPrediction imuPrediction_;
};

}
#endif /* GIF_LANDMARKPREDICTION_HPP_ */
