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

  // Covariances
  const double IrIM_pre_ = 1e-6;
  const double qIM_pre = 1e-6;
  const double MvM_pre = 1e-4;         // IMU acc
  const double MwM_pre = 4e-6;         // IMU gyr
  const double MwM_bias_pre = 1e-8;    // IMU gyr bias
  const double MfM_bias_pre = 1e-8;    // IMU acc bias

  const double IrIM_init = 1e-8;
  const double qIM_init = 1e-2;
  const double MvM_init = 1e-2;
  const double MwM_init = 1e-2;
  const double MwM_bias_init = 1e-2;
  const double MfM_bias_init = 1e-2;

  const double dyn_pos_pre = 1e-1;
  const double dyn_att_pre = 1e-1;
  const double mass_pre = 1e-8;
  const double cT_pre = 1e-8;
  const double cM_pre = 1e-8;
  const double cD_pre = 1e-8;
  const double BrBM_pre = 1e-8;
  const double BrBC_pre = 1e-8;
  const double I_pre = 1e-8;
  const double qMB_pre = 1e-8;

  const double mass_init = 1e-2;
  const double cT_init = 1e-2;
  const double cM_init = 1e-2;
  const double cD_init = 1e-2;
  const double BrBM_init = 1e-2;
  const double BrBC_init = 1e-2;
  const double I_init = 1e-2;
  const double qMB_init = 10;

  const double IrIJ_pre = 1e-8;
  const double qIJ_pre = 1e-8;
  const double MrMC_pre = 1e-8;
  const double qMC_pre = 1e-8;
  const double JrJC_upd = 1e-4;
  const double qJC_upd = 1e-4;

  const double IrIJ_init = 100;
  const double qIJ_init = 10;
  const double MrMC_init = 1;
  const double qMC_init = 10;

  const double zRef_init = 1e6;
  const double zRef_pre = 1e-2;
  const double z_upd = 0.1;

  const bool includeIMU = true;
  const bool includeDyn = true;
  const bool includePose = false;
  const bool includeHeight = true;

  BlindMavFilter():
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
   heightPrediction_(new HeightPrediction("HeightPrediction", {"zRef"}, {"zRef"}))
  {
    include_max_ = !includePose & !includeHeight;

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
    if(includeIMU){
      imuBiasPrediction_->GetNoiseCovarianceBlock("MwM_bias") = GIF::Mat3::Identity()*MwM_bias_pre;
      imuBiasPrediction_->GetNoiseCovarianceBlock("MfM_bias") = GIF::Mat3::Identity()*MfM_bias_pre;
      imu_bias_prediction_id_ = AddResidual(imuBiasPrediction_, GIF::fromSec(10.0), GIF::fromSec(0.0));

      imuaccFindif_->GetNoiseCovarianceBlock("MvM") = GIF::Mat3::Identity()*MvM_pre;
      imuacc_findif_id_ = AddResidual(imuaccFindif_,
                                                      GIF::fromSec(10.0), GIF::fromSec(0.0));

      imurorUpdate_->GetNoiseCovarianceBlock("MwM") = GIF::Mat3::Identity()*MwM_pre;
      imuror_update_id_ = AddResidual(imurorUpdate_,
                                                      GIF::fromSec(10.0), GIF::fromSec(0.0));
    }

    positionFindif_->GetNoiseCovarianceBlock("IrIM") = GIF::Mat3::Identity()*IrIM_pre_;
    position_findif_id_ = AddResidual(positionFindif_,
                                                    GIF::fromSec(10.0), GIF::fromSec(0.0));

    attitudeFindif_->GetNoiseCovarianceBlock("qIM_vec") = GIF::Mat3::Identity()*qIM_pre;
    attitude_findif_id_ = AddResidual(attitudeFindif_,
                                                    GIF::fromSec(10.0), GIF::fromSec(0.0));

    if(includeDyn){
//      mavDynamicResidual_->SetCalibrationFlags(0,1,1,1,1,1,0,1);
      mavDynamicResidual_->SetCalibrationFlags(0,0,0,0,0,0,0,0);

      // Aggressive dataset (identified)
      mavDynamicResidual_->SetMass(1.678);
      mavDynamicResidual_->SetAerodynamicCoefficient(5.69629e-6*1e6,-5.43366e-8*1e6,0.0369976);
      mavDynamicResidual_->SetInertiaDiagonal(1.48542e-2,1.60995e-2,1.34006e-2);
      mavDynamicResidual_->SetComOffset(Vec3(0.000169187, -0.000322526,    0.0324966));
      mavDynamicResidual_->SetImuOffset(Vec3(-0.000736413, 0.00108263,  0.0337457),Quat(1,0,0,0));


      // New datasets (identified)
      mavDynamicResidual_->SetAerodynamicCoefficient(5.90423,-0.0306806,0.0205564);
      mavDynamicResidual_->SetComOffset(Vec3(0.0048764, -6.45603e-05, 0.0116824));
      mavDynamicResidual_->SetImuOffset(Vec3(0.000565769, 0.00741146, 0.0201034),Quat(1,0,0,0));
      // New datasets, better IMU (identified)
      mavDynamicResidual_->SetAerodynamicCoefficient(5.90668,-0.0305815,0.019972);
      mavDynamicResidual_->SetComOffset(Vec3(0.00473559, -0.00016262, 0.0360432));
      mavDynamicResidual_->SetImuOffset(Vec3(0.11471, 0.033082, 0.0239538),Quat(-0.0379607, 0.814306, -0.0509552, 0.576947));


      mavDynamicResidual_->GetNoiseCovarianceBlock("dyn_pos") = GIF::Mat<6>::Identity()*dyn_pos_pre;
      mavDynamicResidual_->GetNoiseCovarianceBlock("dyn_att") = GIF::Mat<6>::Identity()*dyn_att_pre;
      mav_dynamic_residual_id_ = AddResidual(mavDynamicResidual_,
                                             GIF::fromSec(10.0), GIF::fromSec(0.0));

      mavParPrediction_->GetNoiseCovarianceBlock("m") = GIF::Mat<1>::Identity()*mass_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("cT") = GIF::Mat<1>::Identity()*cT_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("cM") = GIF::Mat<1>::Identity()*cM_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("cD") = GIF::Mat<1>::Identity()*cD_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("BrBC") = GIF::Mat<3>::Identity()*BrBC_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("BrBM") = GIF::Mat<3>::Identity()*BrBM_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("I") = GIF::Mat<3>::Identity()*I_pre;
      mavParPrediction_->GetNoiseCovarianceBlock("qMB") = GIF::Mat<3>::Identity()*qMB_pre;
      mavPar_prediction_id_ = AddResidual(mavParPrediction_, GIF::fromSec(10.0), GIF::fromSec(0.0));
    }

    if(includePose){
      poseUpdate_->GetNoiseCovarianceBlock("JrJC") = GIF::Mat<3>::Identity()*JrJC_upd;
      poseUpdate_->GetNoiseCovarianceBlock("qJC") = GIF::Mat<3>::Identity()*qJC_upd;
      pose_IB_update_id_ = AddResidual(poseUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));
      externalCalibPrediction_->GetNoiseCovarianceBlock("IrIJ") = GIF::Mat<3>::Identity()*IrIJ_pre;
      externalCalibPrediction_->GetNoiseCovarianceBlock("qIJ") = GIF::Mat<3>::Identity()*qIJ_pre;
      externalCalibPrediction_->GetNoiseCovarianceBlock("MrMC") = GIF::Mat<3>::Identity()*MrMC_pre;
      externalCalibPrediction_->GetNoiseCovarianceBlock("qMC") = GIF::Mat<3>::Identity()*qMC_pre;
      external_calib_prediction_id_ = AddResidual(externalCalibPrediction_,
                                                  GIF::fromSec(10.0), GIF::fromSec(0.0));
    }
    if(includeHeight){
      heightUpdate_->GetNoiseCovarianceBlock("z") = GIF::Mat<1>::Identity()*z_upd;
      height_update_id_ = AddResidual(heightUpdate_, GIF::fromSec(10.0), GIF::fromSec(0.0));
      heightPrediction_->GetNoiseCovarianceBlock("zRef") = GIF::Mat<1>::Identity()*zRef_pre;
      height_prediction_id_ = AddResidual(heightPrediction_,GIF::fromSec(10.0), GIF::fromSec(0.0));
    }

    std::cout << PrintConnectivity();
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
    bool check = !includeIMU;
    std::shared_ptr<const GIF::ElementVectorBase> accMeas;
    if(includeIMU){
      check = residuals_.at(imuacc_findif_id_).mt_.GetFirst(accMeas);
    }

    if(check){
      startTime_ = t;
      time_ = t;
      inf_.setIdentity();
      GetNoiseInfBlock("IrIM") = GIF::Mat3::Identity()/IrIM_init;
      GetNoiseInfBlock("MvM") = GIF::Mat3::Identity()/MvM_init;
      GetNoiseInfBlock("MwM") = GIF::Mat3::Identity()/MwM_init;
      GetNoiseInfBlock("qIM") = GIF::Mat3::Identity()/qIM_init;
      if(includeIMU){
        GetNoiseInfBlock("MwM_bias") = GIF::Mat3::Identity()/MwM_bias_init;
        GetNoiseInfBlock("MfM_bias") = GIF::Mat3::Identity()/MfM_bias_init;
      }
      if(includeDyn){
        GetNoiseInfBlock("m") = GIF::Mat<1>::Identity()/mass_init;
        GetNoiseInfBlock("cT") = GIF::Mat<1>::Identity()/cT_init;
        GetNoiseInfBlock("cM") = GIF::Mat<1>::Identity()/cM_init;
        GetNoiseInfBlock("cD") = GIF::Mat<1>::Identity()/cD_init;
        GetNoiseInfBlock("BrBC") = GIF::Mat<3>::Identity()/BrBC_init;
        GetNoiseInfBlock("BrBM") = GIF::Mat<3>::Identity()/BrBM_init;
        GetNoiseInfBlock("I") = GIF::Mat<3>::Identity()/I_init;
        GetNoiseInfBlock("qMB") = GIF::Mat<3>::Identity()/qMB_init;
        state_.GetValue<double>("m") = 1.678;
        state_.GetValue<double>("cT") = 5.73311e-6 * 1e6;
        state_.GetValue<double>("cM") = -4.2101e-8 * 1e6;
        state_.GetValue<double>("cD") = 0.018529;
        state_.GetValue<GIF::Vec3>("BrBC") = GIF::Vec3(0.0241, -0.0140, -0.0372);
        state_.GetValue<GIF::Vec3>("BrBM") = GIF::Vec3(-0.0275, 0.0158, 0.0377);
        state_.GetValue<GIF::Vec3>("I") = GIF::Vec3(1.48542e-2,1.60995e-2,1.34006e-2);
        state_.GetValue<GIF::Quat>("qMB") = GIF::Quat(1,0,0,0);
        state_.GetValue<GIF::Quat>("qMB") = GIF::Quat(-0.00464142246184, 0.813813125739, 0.00377470375119, 0.581095865743);

        state_.GetValue<double>("cT") = 5.90387;
  //      state_.GetValue<double>("cM") = -3.0513e-8 * 1e6;
  //      state_.GetValue<double>("cD") = 0.018609;
  //      state_.GetValue<GIF::Vec3>("BrBC") = GIF::Vec3(0.00489592, -4.94558e-05, 0.011715);
  //      state_.GetValue<GIF::Vec3>("BrBM") = GIF::Vec3(-0.00109527, 0.00606688, 0.0202162);
      }
      if(includePose){
        GetNoiseInfBlock("IrIJ") = GIF::Mat<3>::Identity()/IrIJ_init;
        GetNoiseInfBlock("qIJ") = GIF::Mat<3>::Identity()/qIJ_init;
        GetNoiseInfBlock("MrMC") = GIF::Mat<3>::Identity()/MrMC_init;
        GetNoiseInfBlock("qMC") = GIF::Mat<3>::Identity()/qMC_init;
        state_.GetValue<GIF::Vec3>("IrIJ").setZero();
        state_.GetValue<GIF::Quat>("qIJ") = GIF::Quat(-0.155377, 0.000959049, 0.00422654, 0.987846);
//        state_.GetValue<GIF::Quat>("qIJ") = GIF::Quat(1,0,0,0);
        state_.GetValue<GIF::Vec3>("MrMC").setZero();
        state_.GetValue<GIF::Quat>("qMC") = GIF::Quat(-0.00464142246184, 0.813813125739, 0.00377470375119, 0.581095865743);
//        state_.GetValue<GIF::Quat>("qMC") = GIF::Quat(1,0,0,0);
      }
      if(includeHeight){
        GetNoiseInfBlock("zRef") = GIF::Mat<1>::Identity()/zRef_init;
        state_.GetValue<double>("zRef") = 0;
      }

      if(includeIMU){
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
};

}
#endif /* GIF_BLINDMAVFILTER_HPP_ */
