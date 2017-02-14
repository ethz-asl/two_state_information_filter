#ifndef GIF_BLINDMAVFILTER_HPP_
#define GIF_BLINDMAVFILTER_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/filter.h"
#include "generalized_information_filter/residuals/height-update.h"
#include "generalized_information_filter/residuals/mav-dynamic-findif.h"
#include "generalized_information_filter/residuals/pose-update.h"
#include "generalized_information_filter/residuals/random-walk-prediction.h"
#include "generalized_information_filter/residuals/robocentric/attitude-findif.h"
#include "generalized_information_filter/residuals/robocentric/imuacc-findif.h"
#include "generalized_information_filter/residuals/robocentric/imuror-update.h"
#include "generalized_information_filter/residuals/robocentric/position-findif.h"

namespace GIF {

/*! \brief Filter for fusing proprioception on an MAV
 *         It is meant for fusing the following sensor modalities:
 *         - IMU
 *         - Rotor speeds
 *         - Pressure sensor
 *         It provides some additional calibration capabilities.
 */
class BlindMavFilter: public GIF::Filter{
 public:
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<GIF::Vec3,GIF::Vec3>> ImuBiasPrediction;
  typedef GIF::ImuaccFindif ImuaccFindif;
  typedef GIF::PositionFindif PositionFindif;
  typedef GIF::AttitudeFindif AttitudeFindif;
  typedef GIF::ImurorUpdate ImurorUpdate;
  typedef GIF::MavDynamicFindif<6> MavDynamicResidual;
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<double,double,double,double,GIF::Vec3,
                                                     GIF::Vec3,GIF::Vec3,GIF::Quat>> MavParPrediction;
  typedef GIF::PoseUpdate<true,true> PoseUpdate;
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<GIF::Vec3,GIF::Quat,
                                                     GIF::Vec3,GIF::Quat>> ExternalCalibPrediction;
  typedef GIF::HeightUpdate HeightUpdate;
  typedef GIF::RandomWalkPrediction<GIF::ElementPack<double>> HeightPrediction;
  std::shared_ptr<ImuBiasPrediction> imuBiasPrediction_;
  std::shared_ptr<ImuaccFindif> imuaccFindif_;
  std::shared_ptr<PositionFindif> positionFindif_;
  std::shared_ptr<AttitudeFindif> attitudeFindif_;
  std::shared_ptr<ImurorUpdate> imurorUpdate_;
  std::shared_ptr<MavDynamicResidual> mavDynamicResidual_;
  std::shared_ptr<MavParPrediction> mavParPrediction_;
  std::shared_ptr<PoseUpdate> poseUpdate_;
  std::shared_ptr<ExternalCalibPrediction> externalCalibPrediction_;
  std::shared_ptr<HeightUpdate> heightUpdate_;
  std::shared_ptr<HeightPrediction> heightPrediction_;
  int imu_bias_prediction_id_;
  int imuacc_findif_id_;
  int position_findif_id_;
  int attitude_findif_id_;
  int imuror_update_id_;
  int mav_dynamic_residual_id_;
  int mavPar_prediction_id_;
  int pose_IB_update_id_;
  int external_calib_prediction_id_;
  int height_update_id_;
  int height_prediction_id_;

  // Init Covariances
  double MwM_bias_init_;
  double MfM_bias_init_;
  double IrIM_init_;
  double qIM_init_;
  double MvM_init_;
  double MwM_init_;

  double mass_init_;
  double I_init_;
  double cT_init_;
  double cM_init_;
  double cD_init_;
  double BrBM_init_;
  double BrBC_init_;
  double qMB_init_;

  double IrIJ_init_;
  double qIJ_init_;
  double MrMC_init_;
  double qMC_init_;
  double zRef_init_;

  const bool useImu_;
  const bool useDyn_;
  const bool usePose_;
  const bool useHeight_;

  // Init State
  ElementVector initState_;

  BlindMavFilter(bool useImu, bool useDyn, bool useHeight, bool usePose):
    imuBiasPrediction_(new ImuBiasPrediction("ImuBiasPrediction", {"MwM_bias", "MfM_bias"},
                                             {"MwM_bias", "MfM_bias"})),
    imuaccFindif_(new ImuaccFindif("ImuaccFindif", {"MvM"}, {"MvM", "MwM", "MfM_bias", "qIM"},
                                   {"MvM"}, {"MvM"})),
    positionFindif_(new PositionFindif("PositionFindif", {"IrIM"}, {"IrIM", "MvM", "qIM"},
                                       {"IrIM"}, {"IrIM"})),
    attitudeFindif_(new AttitudeFindif("AttitudeFindif", {"qIM_vec"}, {"qIM", "MwM"}, {"qIM"},
                                       {"qIM_vec"})),
    imurorUpdate_(new ImurorUpdate("ImurorUpdate", {"MwM"}, {}, {"MwM", "MwM_bias"}, {"MwM"})),
    mavDynamicResidual_(new MavDynamicResidual("MavDynamicResidual", {"dyn_pos", "dyn_att"},
        {"MvM", "MwM", "qIM", "m", "cT", "cM", "cD", "BrBM", "BrBC", "I", "qMB"},
        {"MvM", "MwM"}, {"dyn_pos", "dyn_att"})),
    mavParPrediction_(new MavParPrediction("MavParPrediction",
                                           {"m", "cT", "cM", "cD", "BrBM", "BrBC", "I", "qMB"},
                                           {"m", "cT", "cM", "cD", "BrBM", "BrBC", "I", "qMB"})),
    poseUpdate_(new PoseUpdate("PoseUpdate", {"JrJC", "qJC"},
                               {"IrIM", "qIM", "IrIJ", "qIJ", "MrMC", "qMC"}, {"JrJC", "qJC"})),
    externalCalibPrediction_(new ExternalCalibPrediction("ExternalCalibPrediction",
                                                         {"IrIJ", "qIJ", "MrMC", "qMC"},
                                                         {"IrIJ", "qIJ", "MrMC", "qMC"})),
   heightUpdate_(new HeightUpdate("HeightUpdate", {"z"}, {"zRef", "IrIM"}, {"z"})),
   heightPrediction_(new HeightPrediction("HeightPrediction", {"zRef"}, {"zRef"})),
   initState_(stateDefinition_),
   useImu_(useImu), useDyn_(useDyn), useHeight_(useHeight), usePose_(usePose)
  {
    include_max_ = !usePose_ & !useHeight_;

    imu_bias_prediction_id_ = -1;
    imuacc_findif_id_ = -1;
    position_findif_id_ = -1;
    attitude_findif_id_ = -1;
    imuror_update_id_ = -1;
    mav_dynamic_residual_id_ = -1;
    mavPar_prediction_id_ = -1;
    pose_IB_update_id_ = -1;
    external_calib_prediction_id_ = -1;
    height_update_id_ = -1;
    height_prediction_id_ = -1;

    // Imu Bias Prediction Properties
    if(useImu_){
      imu_bias_prediction_id_ = AddResidual(imuBiasPrediction_, GIF::fromSec(10.0),
                                            GIF::fromSec(0.0));
      imuacc_findif_id_ = AddResidual(imuaccFindif_, GIF::fromSec(10.0), GIF::fromSec(0.0));
      imuror_update_id_ = AddResidual(imurorUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));
    }
    position_findif_id_ = AddResidual(positionFindif_, GIF::fromSec(10.0), GIF::fromSec(0.0));
    attitude_findif_id_ = AddResidual(attitudeFindif_, GIF::fromSec(10.0), GIF::fromSec(0.0));

    if(useDyn_){
      mav_dynamic_residual_id_ = AddResidual(mavDynamicResidual_, GIF::fromSec(10.0),
                                             GIF::fromSec(0.0));
      mavPar_prediction_id_ = AddResidual(mavParPrediction_, GIF::fromSec(10.0), GIF::fromSec(0.0));
    }

    if(usePose_){
      pose_IB_update_id_ = AddResidual(poseUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));
      external_calib_prediction_id_ = AddResidual(externalCalibPrediction_,GIF::fromSec(10.0),
                                                  GIF::fromSec(0.0));
    }
    if(useHeight_){
      height_update_id_ = AddResidual(heightUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));
      height_prediction_id_ = AddResidual(heightPrediction_,GIF::fromSec(10.0), GIF::fromSec(0.0));
    }
    std::cout << PrintConnectivity();
    initState_.Construct();

    SetIterationParameters(10,0.1);
    SetDynCalibrationFlags(0,0,0,0,0,0,0,0);

    // Covariance parameter
    SetImuCovParameter(1e-4, 4e-6, 1e-8, 1e-8, 1e-2, 1e-2);
    SetIntegrationCovParameter(1e-6,1e-6, 1e-8, 1e-2, 1e-2, 1e-2);
    SetDynamicCovParameter(1e-1,1e-1,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,
                           1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,10);
    SetPoseCovParameter(1e-8,1e-8,1e-8,1e-8,1e-4,1e-4,100,10,1,10);
    SetHeightCovParameter(1e-2,0.1,1e6);

    // Dynamic parameter (Firefly - bluebird)
    SetMass(1.678);
    SetAerodynamicCoefficient(5.73311,-0.042101,0.018529);
    SetInertiaDiagonal(1.48542e-2,1.60995e-2,1.34006e-2);
    SetDynComOffset(Vec3(0.0241, -0.0140, -0.0372));
    SetDynImuOffset(Vec3(-0.0275, 0.0158, 0.0377),Quat(1,0,0,0));

    // Pose parameter
    SetInertialAlignment(Vec3(0,0,0),Quat(1,0,0,0));
    SetBodyAlignment(Vec3(0,0,0),Quat(1,0,0,0));
  }
  virtual ~BlindMavFilter(){};
  void PreProcess(){
  };
  void PostProcess(){
  };
  virtual void ComputeLinearizationPoint(const GIF::TimePoint& t){
    GIF::ElementVectorBase::CPtr accMeas;
    GIF::ElementVectorBase::CPtr gyrMeas;

    curLinState_ = state_;
    if(residuals_.at(imuacc_findif_id_).mt_.GetMeasurement(t, accMeas)
        && residuals_.at(imuror_update_id_).mt_.GetMeasurement(t, gyrMeas)){
      std::shared_ptr<const GIF::AccMeas> acc = std::dynamic_pointer_cast<const GIF::AccMeas>(accMeas);
      std::shared_ptr<const GIF::RorMeas> gyr = std::dynamic_pointer_cast<const GIF::RorMeas>(gyrMeas);

      const double dt = GIF::toSec(t-time_);
      curLinState_.GetValue<GIF::Vec3>("MwM") = gyr->MwM_ - state_.GetValue<GIF::Vec3>("MwM_bias");
      curLinState_.GetValue<GIF::Vec3>("IrIM") = state_.GetValue<GIF::Vec3>("IrIM")
          + dt * state_.GetValue<GIF::Quat>("qIM").rotate(state_.GetValue<GIF::Vec3>("MvM"));
      curLinState_.GetValue<GIF::Vec3>("MvM") =
          (GIF::Mat3::Identity() - GIF::gSM(dt * state_.GetValue<GIF::Vec3>("MwM"))) * state_.GetValue<GIF::Vec3>("MvM")
          + dt * (acc->MfM_ - state_.GetValue<GIF::Vec3>("MfM_bias")
                  + state_.GetValue<GIF::Quat>("qIM").inverseRotate(GIF::Vec3(0,0,-9.81)));
      GIF::Quat dQ = dQ.exponentialMap(dt * state_.GetValue<GIF::Vec3>("MwM"));
      curLinState_.GetValue<GIF::Quat>("qIM") = state_.GetValue<GIF::Quat>("qIM") * dQ;
    }
  }
  void Init(const GIF::TimePoint& t = GIF::TimePoint::min(), const GIF::ElementVectorBase* initState = nullptr) {
    bool check = !useImu_;
    std::shared_ptr<const GIF::ElementVectorBase> accMeas;
    if(useImu_){
      check = residuals_.at(imuacc_findif_id_).mt_.GetFirst(accMeas);
    }

    if(check){
      startTime_ = t;
      time_ = t;
      inf_.setIdentity();
      state_ = initState_;
      GetNoiseInfBlock("IrIM") = GIF::Mat3::Identity()/IrIM_init_;
      GetNoiseInfBlock("qIM") = GIF::Mat3::Identity()/qIM_init_;
      GetNoiseInfBlock("MvM") = GIF::Mat3::Identity()/MvM_init_;
      GetNoiseInfBlock("MwM") = GIF::Mat3::Identity()/MwM_init_;
      if(useImu_){
        GetNoiseInfBlock("MwM_bias") = GIF::Mat3::Identity()/MwM_bias_init_;
        GetNoiseInfBlock("MfM_bias") = GIF::Mat3::Identity()/MfM_bias_init_;
      }
      if(useDyn_){
        GetNoiseInfBlock("m") = GIF::Mat<1>::Identity()/mass_init_;
        GetNoiseInfBlock("cT") = GIF::Mat<1>::Identity()/cT_init_;
        GetNoiseInfBlock("cM") = GIF::Mat<1>::Identity()/cM_init_;
        GetNoiseInfBlock("cD") = GIF::Mat<1>::Identity()/cD_init_;
        GetNoiseInfBlock("BrBC") = GIF::Mat<3>::Identity()/BrBC_init_;
        GetNoiseInfBlock("BrBM") = GIF::Mat<3>::Identity()/BrBM_init_;
        GetNoiseInfBlock("I") = GIF::Mat<3>::Identity()/I_init_;
        GetNoiseInfBlock("qMB") = GIF::Mat<3>::Identity()/qMB_init_;
      }
      if(usePose_){
        GetNoiseInfBlock("IrIJ") = GIF::Mat<3>::Identity()/IrIJ_init_;
        GetNoiseInfBlock("qIJ") = GIF::Mat<3>::Identity()/qIJ_init_;
        GetNoiseInfBlock("MrMC") = GIF::Mat<3>::Identity()/MrMC_init_;
        GetNoiseInfBlock("qMC") = GIF::Mat<3>::Identity()/qMC_init_;
      }
      if(useHeight_){
        GetNoiseInfBlock("zRef") = GIF::Mat<1>::Identity()/zRef_init_;
      }

      if(useImu_){
        // Use accelerometer to estimate initial attitude
        GIF::Vec3 unitZ(0,0,1);
        const GIF::Vec3 MfM = std::dynamic_pointer_cast<const GIF::AccMeas>(accMeas)->MfM_;
        // Align with gravity
        if(MfM.norm()>1e-6){
          GIF::Quat q;
          q.setFromVectors(state_.GetValue<GIF::Quat>("qIM").rotate(MfM),unitZ);
          state_.GetValue<GIF::Quat>("qIM") = q*state_.GetValue<GIF::Quat>("qIM");
        }
      }
      // TODO: use first height measurement for initializing zRef
      LOG(INFO) << "Initializing state at t = " << GIF::Print(t) << std::endl;
      is_initialized_ = true;
    }
  }
  void SetImuCovParameter(double MvM_pre, double MwM_pre, double MwM_bias_pre, double MfM_bias_pre,
                       double MwM_bias_init, double MfM_bias_init){
    imuaccFindif_->GetNoiseCovarianceBlock("MvM") = GIF::Mat3::Identity()*MvM_pre;
    imurorUpdate_->GetNoiseCovarianceBlock("MwM") = GIF::Mat3::Identity()*MwM_pre;
    imuBiasPrediction_->GetNoiseCovarianceBlock("MwM_bias") = GIF::Mat3::Identity()*MwM_bias_pre;
    imuBiasPrediction_->GetNoiseCovarianceBlock("MfM_bias") = GIF::Mat3::Identity()*MfM_bias_pre;
    MwM_bias_init_ = MwM_bias_init;
    MfM_bias_init_ = MfM_bias_init;
  }
  void SetIntegrationCovParameter(double IrIM_pre, double qIM_pre,
                               double IrIM_init, double qIM_init, double MvM_init, double MwM_init){
    positionFindif_->GetNoiseCovarianceBlock("IrIM") = GIF::Mat3::Identity()*IrIM_pre;
    attitudeFindif_->GetNoiseCovarianceBlock("qIM_vec") = GIF::Mat3::Identity()*qIM_pre;
    IrIM_init_ = IrIM_init;
    qIM_init_ = qIM_init;
    MvM_init_ = MvM_init;
    MwM_init_ = MwM_init;
  }
  void SetDynamicCovParameter(double dyn_pos_pre, double dyn_att_pre, double mass_pre, double I_pre,
                           double cT_pre, double cM_pre, double cD_pre, double BrBC_pre,
                           double BrBM_pre, double qMB_pre,
                           double mass_init, double I_init, double cT_init, double cM_init,
                           double cD_init, double BrBM_init, double BrBC_init, double qMB_init){
    mavDynamicResidual_->GetNoiseCovarianceBlock("dyn_pos") = GIF::Mat<6>::Identity()*dyn_pos_pre;
    mavDynamicResidual_->GetNoiseCovarianceBlock("dyn_att") = GIF::Mat<6>::Identity()*dyn_att_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("m") = GIF::Mat<1>::Identity()*mass_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("cT") = GIF::Mat<1>::Identity()*cT_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("cM") = GIF::Mat<1>::Identity()*cM_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("cD") = GIF::Mat<1>::Identity()*cD_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("BrBC") = GIF::Mat<3>::Identity()*BrBC_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("BrBM") = GIF::Mat<3>::Identity()*BrBM_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("I") = GIF::Mat<3>::Identity()*I_pre;
    mavParPrediction_->GetNoiseCovarianceBlock("qMB") = GIF::Mat<3>::Identity()*qMB_pre;

    mass_init_ = mass_init;
    I_init_ = I_init;
    cT_init_ = cT_init;
    cM_init_ = cM_init;
    cD_init_ = cD_init;
    BrBM_init_ = BrBM_init;
    BrBC_init_ = BrBC_init;
    qMB_init_ = qMB_init;
  }
  void SetPoseCovParameter(double IrIJ_pre, double qIJ_pre, double MrMC_pre, double qMC_pre,
                        double JrJC_upd, double qJC_upd,
                        double IrIJ_init, double qIJ_init, double MrMC_init, double qMC_init){

    poseUpdate_->GetNoiseCovarianceBlock("JrJC") = GIF::Mat<3>::Identity()*JrJC_upd;
    poseUpdate_->GetNoiseCovarianceBlock("qJC") = GIF::Mat<3>::Identity()*qJC_upd;
    externalCalibPrediction_->GetNoiseCovarianceBlock("IrIJ") = GIF::Mat<3>::Identity()*IrIJ_pre;
    externalCalibPrediction_->GetNoiseCovarianceBlock("qIJ") = GIF::Mat<3>::Identity()*qIJ_pre;
    externalCalibPrediction_->GetNoiseCovarianceBlock("MrMC") = GIF::Mat<3>::Identity()*MrMC_pre;
    externalCalibPrediction_->GetNoiseCovarianceBlock("qMC") = GIF::Mat<3>::Identity()*qMC_pre;

    IrIJ_init_ = IrIJ_init;
    qIJ_init_ = qIJ_init;
    MrMC_init_ = MrMC_init;
    qMC_init_ = qMC_init;
  }
  void SetHeightCovParameter(double zRef_pre, double z_upd, double zRef_init){
    heightUpdate_->GetNoiseCovarianceBlock("z") = GIF::Mat<1>::Identity()*z_upd;
    heightPrediction_->GetNoiseCovarianceBlock("zRef") = GIF::Mat<1>::Identity()*zRef_pre;

    zRef_init_ = zRef_init;
  }
  void SetDynImuOffset(Vec3 BrBM, Quat qMB){
    mavDynamicResidual_->SetImuOffset(BrBM,qMB);
    if(useDyn_){
      initState_.GetValue<GIF::Vec3>("BrBM") = BrBM;
      initState_.GetValue<GIF::Quat>("qMB") = qMB;
    }
  }
  void SetDynComOffset(Vec3 BrBC){
    mavDynamicResidual_->SetComOffset(BrBC);
    if(useDyn_){
      initState_.GetValue<GIF::Vec3>("BrBC") = BrBC;
    }
  }
  void SetMass(double m){
    mavDynamicResidual_->SetMass(m);
    if(useDyn_){
      initState_.GetValue<double>("m") = m;
    }
  }
  void SetAerodynamicCoefficient(double cT, double cM, double cD){
    mavDynamicResidual_->SetAerodynamicCoefficient(cT,cM,cD);
    if(useDyn_){
      initState_.GetValue<double>("cT") = cT;
      initState_.GetValue<double>("cM") = cM;
      initState_.GetValue<double>("cD") = cD;
    }
  }
  void SetInertiaDiagonal(double Ix, double Iy, double Iz){
    mavDynamicResidual_->SetInertiaDiagonal(Ix,Iy,Iz);
    if(useDyn_){
      initState_.GetValue<GIF::Vec3>("I") = GIF::Vec3(Ix,Iy,Iz);
    }
  }
  void SetDynCalibrationFlags(bool m, bool cT, bool cM, bool cD, bool BrBM, bool BrBC, bool I, bool qMB){
    mavDynamicResidual_->SetCalibrationFlags(m, cT, cM, cD, BrBM, BrBC, I, qMB);
  }
  void SetInertialAlignment(const Vec3& IrIJ, const Quat& qIJ){
    poseUpdate_->SetInertialAlignment(IrIJ,qIJ);
    if(usePose_){
      initState_.GetValue<GIF::Vec3>("IrIJ") = IrIJ;
      initState_.GetValue<GIF::Quat>("qIJ") = qIJ;
    }
  }
  void SetBodyAlignment(const Vec3& BrBC, const Quat& qBC){
    poseUpdate_->SetInertialAlignment(BrBC,qBC);
    if(usePose_){
      initState_.GetValue<GIF::Vec3>("MrMC") = BrBC;
      initState_.GetValue<GIF::Quat>("qMC") = qBC;
    }
  }
  void SetIterationParameters(int num_iter, double iter_th){
    num_iter_ = num_iter;
    iter_th_ = iter_th;
  }
};

}
#endif /* GIF_BLINDMAVFILTER_HPP_ */
