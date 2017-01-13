#ifndef GIF_IMUPOSEFILTER_HPP_
#define GIF_IMUPOSEFILTER_HPP_

#include "generalized_information_filter/residuals/robocentric/attitude-findif.h"
#include "generalized_information_filter/residuals/robocentric/imuacc-findif.h"
#include "generalized_information_filter/residuals/robocentric/imuror-update.h"
#include "generalized_information_filter/residuals/robocentric/position-findif.h"
#include "generalized_information_filter/common.h"
#include "generalized_information_filter/residuals/random-walk-prediction.h"
#include "generalized_information_filter/residuals/pose-update.h"
#include "generalized_information_filter/filter.h"

namespace GIF {

/*! \brief Filter for fusing IMU measurements with pose measurements
 *         An arbitrary number of IMUs can be included
 */
template<int NumImu, bool doInertialAlignment, bool doBodyAlignment>
class ImuPoseFilter: public GIF::Filter{
 public:
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<GIF::Vec3,GIF::Vec3>> ImuBiasPrediction;
  typedef GIF::ImuaccFindif ImuaccFindif;
  typedef GIF::PositionFindif PositionFindif;
  typedef GIF::AttitudeFindif AttitudeFindif;
  typedef GIF::ImurorUpdate ImurorUpdate;
  typedef GIF::PoseUpdate<doInertialAlignment,doBodyAlignment> PoseUpdate;
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<GIF::Vec3,GIF::Quat>> PosePrediction;
  std::shared_ptr<ImuBiasPrediction> imuBiasPrediction_[NumImu];
  std::shared_ptr<ImuaccFindif> imuaccFindif_[NumImu];
  std::shared_ptr<PositionFindif> positionFindif_;
  std::shared_ptr<AttitudeFindif> attitudeFindif_;
  std::shared_ptr<ImurorUpdate> imurorUpdate_[NumImu];
  std::shared_ptr<PoseUpdate> poseUpdate_;
  std::shared_ptr<PosePrediction> inertialCalibPrediction_;
  std::shared_ptr<PosePrediction> bodyCalibPrediction_;
  int imu_bias_prediction_id_[NumImu];
  int robcen_imuacc_findif_residual_id_[NumImu];
  int robcen_velocity_findif_residual_id_;
  int robcen_rotrate_findif_residual_id_;
  int robcen_imuror_update_id_[NumImu];
  int pose_update_id_;
  int inertial_calib_prediction_id_;
  int body_calib_prediction_id_;

  // Covariances
  const double IrIM_pre = 1e-6;
  const double qIM_pre = 1e-6;
  const double MvM_pre = 1e-2;         // IMU acc
  const double MwM_bias_pre = 1e-8;    // IMU gyr bias
  const double MfM_bias_pre = 1e-8;    // IMU acc bias
  const double MwM_pre = 4e-6;         // IMU gyr

  const double IrIJ_pre = 1e-8;
  const double qIJ_pre = 1e-8;
  const double MrMC_pre = 1e-8;
  const double qMC_pre = 1e-8;
  const double JrJC_upd = 1e-4;
  const double qJC_upd = 1e-4;

  const double IrIM_init = 1e-8;
  const double MvM_init = 1e-2;
  const double MwM_init = 1e-2;
  const double MwM_bias_init = 1e-2;
  const double MfM_bias_init = 1e-2;
  const double qIM_init = 1e-2;

  const double IrIJ_init = 100;
  const double qIJ_init = 10;
  const double MrMC_init = 1;
  const double qMC_init = 1;

  const double MvM_huber_ = 0.2;
  const double MwM_huber_ = 0.2;
  const double JrJC_huber_ = 0.1;

  ImuPoseFilter(){
    positionFindif_.reset(new PositionFindif(
        "PositionFindif", {"IrIM"}, {"IrIM", "MvM", "qIM"}, {"IrIM"}, {"IrIM"}));
    positionFindif_->GetNoiseCovarianceBlock("IrIM") = GIF::Mat3::Identity()*IrIM_pre;
    robcen_velocity_findif_residual_id_ = AddResidual(positionFindif_,
                                                      GIF::fromSec(10.0), GIF::fromSec(0.0));

    attitudeFindif_.reset(new AttitudeFindif(
        "AttitudeFindif", {"qIM_vec"}, {"qIM", "MwM"}, {"qIM"}, {"qIM_vec"}));
    attitudeFindif_->GetNoiseCovarianceBlock("qIM_vec") = GIF::Mat3::Identity()*qIM_pre;
    robcen_rotrate_findif_residual_id_ = AddResidual(attitudeFindif_,
                                                     GIF::fromSec(10.0), GIF::fromSec(0.0));

    inertial_calib_prediction_id_ = -1;
    body_calib_prediction_id_ = -1;
    std::array<std::string,2+2*doBodyAlignment+2*doInertialAlignment> poseStateName;
    poseStateName[0] = "IrIM";
    poseStateName[1] = "qIM";
    if(doInertialAlignment){
      poseStateName[0+2*doInertialAlignment] = "IrIJ";
      poseStateName[1+2*doInertialAlignment] = "qIJ";
      inertialCalibPrediction_.reset(new PosePrediction("InertialCalibPrediction",
        {"IrIJ", "qIJ"}, {"IrIJ", "qIJ"}));
      inertialCalibPrediction_->GetNoiseCovarianceBlock("IrIJ") = GIF::Mat<3>::Identity()*IrIJ_pre;
      inertialCalibPrediction_->GetNoiseCovarianceBlock("qIJ") = GIF::Mat<3>::Identity()*qIJ_pre;
      inertial_calib_prediction_id_ = AddResidual(inertialCalibPrediction_,
                                                  GIF::fromSec(10.0), GIF::fromSec(0.0));
    }
    if(doBodyAlignment){
      poseStateName[0+2*doInertialAlignment+2*doBodyAlignment] = "MrMC";
      poseStateName[1+2*doInertialAlignment+2*doBodyAlignment] = "qMC";
      bodyCalibPrediction_.reset(new PosePrediction("BodyCalibPrediction",
                                                    {"MrMC", "qMC"}, {"MrMC", "qMC"}));
      bodyCalibPrediction_->GetNoiseCovarianceBlock("MrMC") = GIF::Mat<3>::Identity()*MrMC_pre;
      bodyCalibPrediction_->GetNoiseCovarianceBlock("qMC") = GIF::Mat<3>::Identity()*qMC_pre;
      body_calib_prediction_id_ = AddResidual(bodyCalibPrediction_,
                                              GIF::fromSec(10.0), GIF::fromSec(0.0));
    }
    poseUpdate_.reset(new PoseUpdate(
        "PoseUpdate", {"JrJC", "qJC"}, poseStateName, {"JrJC", "qJC"}));
    poseUpdate_->GetNoiseCovarianceBlock("JrJC") = GIF::Mat<3>::Identity()*JrJC_upd;
    poseUpdate_->GetNoiseCovarianceBlock("qJC") = GIF::Mat<3>::Identity()*qJC_upd;
    poseUpdate_->SetBodyAlignment(GIF::Vec3(0.0204726667725, -0.0346315780816, 0.0284436991268),GIF::Quat(-0.999963138703, 0.00313250259206, 0.00788912255922, -0.00129244276308));
    poseUpdate_->SetAttitudeFlag(false);
    poseUpdate_->SetHuberTh(JrJC_huber_);
    pose_update_id_ = AddResidual(poseUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));

    for(int i=0;i<NumImu;i++){
      std::string num = std::to_string(i);

      imuBiasPrediction_[i].reset(new ImuBiasPrediction("ImuBiasPrediction" + num,
          {"MwM_bias" + num, "MfM_bias" + num}, {"MwM_bias" + num, "MfM_bias" + num}));
      imuBiasPrediction_[i]->GetNoiseCovarianceBlock("MwM_bias" + num) = GIF::Mat3::Identity()*MwM_bias_pre;
      imuBiasPrediction_[i]->GetNoiseCovarianceBlock("MfM_bias" + num) = GIF::Mat3::Identity()*MfM_bias_pre;
      imu_bias_prediction_id_[i] = AddResidual(imuBiasPrediction_[i], GIF::fromSec(10.0), GIF::fromSec(0.0));

      imuaccFindif_[i].reset(new ImuaccFindif(
          "ImuaccFindif" + num, {"MvM" + num},
          {"MvM", "MwM", "MfM_bias" + num, "qIM"}, {"MvM"}, {"MvM" + num}));
      imuaccFindif_[i]->GetNoiseCovarianceBlock("MvM" + num) = GIF::Mat3::Identity()*MvM_pre;
      imuaccFindif_[i]->SetHuberTh(MvM_huber_);
      robcen_imuacc_findif_residual_id_[i] = AddResidual(imuaccFindif_[i],
                                                         GIF::fromSec(10.0), GIF::fromSec(0.0));

      imurorUpdate_[i].reset(new ImurorUpdate("ImurorUpdate" + num,
          {"MwM" + num}, {}, {"MwM", "MwM_bias" + num}, {"MwM" + num}));
      imurorUpdate_[i]->GetNoiseCovarianceBlock("MwM" + num) = GIF::Mat3::Identity()*MwM_pre;
      imurorUpdate_[i]->SetHuberTh(MwM_huber_);
      robcen_imuror_update_id_[i] = AddResidual(imurorUpdate_[i],
                                                GIF::fromSec(10.0), GIF::fromSec(0.0));
    }

    std::cout << PrintConnectivity();
  }
  virtual ~ImuPoseFilter(){};
  void Init(const GIF::TimePoint& t = GIF::TimePoint::min(), const GIF::ElementVectorBase* initState = nullptr) {
    startTime_ = t;
    time_ = t;
    inf_.setIdentity();
    GetNoiseInfBlock("IrIM") = GIF::Mat3::Identity()/IrIM_init;
    GetNoiseInfBlock("MvM") = GIF::Mat3::Identity()/MvM_init;
    GetNoiseInfBlock("MwM") = GIF::Mat3::Identity()/MwM_init;
    GetNoiseInfBlock("qIM") = GIF::Mat3::Identity()/qIM_init;
    for(int i=0;i<NumImu;i++){
      std::string num = std::to_string(i);
      GetNoiseInfBlock("MwM_bias" + num) = GIF::Mat3::Identity()/MwM_bias_init;
      GetNoiseInfBlock("MfM_bias" + num) = GIF::Mat3::Identity()/MfM_bias_init;
    }
    if(doInertialAlignment){
      GetNoiseInfBlock("IrIJ") = GIF::Mat<3>::Identity()/IrIJ_init;
      GetNoiseInfBlock("qIJ") = GIF::Mat<3>::Identity()/qIJ_init;
    }
    if(doBodyAlignment){
      GetNoiseInfBlock("MrMC") = GIF::Mat<3>::Identity()/MrMC_init;
      GetNoiseInfBlock("qMC") = GIF::Mat<3>::Identity()/qMC_init;
    }

    // Use accelerometer to estimate initial attitude
    GIF::Vec3 unitZ(0,0,1);
    std::shared_ptr<const GIF::ElementVectorBase> accMeas;
    std::shared_ptr<const GIF::ElementVectorBase> poseMeas;
    if(residuals_.at(robcen_imuacc_findif_residual_id_[0]).mt_.GetFirst(accMeas)){
      const GIF::Vec3 MfM = std::dynamic_pointer_cast<const GIF::AccMeas>(accMeas)->MfM_;
      // Align with gravity
      if(MfM.norm()>1e-6){
        GIF::Quat q;
        q.setFromVectors(state_.GetValue<GIF::Quat>("qIM").rotate(MfM),unitZ);
        state_.GetValue<GIF::Quat>("qIM") = q.inverted()*state_.GetValue<GIF::Quat>("qIM");
      }
    }

    LOG(INFO) << "Initializing state at t = " << GIF::Print(t) << std::endl;
    is_initialized_ = true;
  }
};

}
#endif /* GIF_IMUPOSEFILTER_HPP_ */
