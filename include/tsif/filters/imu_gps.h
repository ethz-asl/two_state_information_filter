#ifndef TSIF_IMU_GPS_H_
#define TSIF_IMU_GPS_H_

#include "tsif/residuals/accelerometer_prediction.h"
#include "tsif/residuals/attitude_findif.h"
#include "tsif/residuals/gyroscope_update.h"
#include "tsif/residuals/pose_update.h"
#include "tsif/residuals/position_findif.h"
#include "tsif/residuals/random_walk.h"
#include "tsif/utils/common.h"
#include "tsif/filter.h"

namespace tsif{

class ImuGpsFilter: public Filter<PositionFindif<0,0,2,1>,
                                  AttitudeFindif<0,1,3>,
                                  AccelerometerPrediction<0,2,1,3,4>,
                                  GyroscopeUpdate<0,3,5>,
                                  RandomWalk<Element<Vec3,4>>,
                                  RandomWalk<Element<Vec3,5>>,
                                  PositionUpdate<0,0,1,-6,-7,-8>,
                                  RandomWalk<Element<Vec3,-6>,Element<Quat,-7>,Element<Vec3,-8>>,
                                  AttitudeUpdate<-1,1,-7,-9>>{
 public:
  ImuGpsFilter(){
    max_iter_ = 10;
    th_iter_ = 0.1;
    std::get<0>(residuals_).w_ = 1e4; // PositionFindif
    std::get<1>(residuals_).w_ = 1e4; // AttitudeFindif
    std::get<2>(residuals_).w_ = 1e2; // AccelerometerPrediction
    std::get<3>(residuals_).w_ = 1/(2e-3); // GyroscopeUpdate
    std::get<4>(residuals_).w_ = 1e4; // AccelerometerBiasPrediction
    std::get<5>(residuals_).w_ = 1e4; // GyroscopeBiasPrediction
    std::get<6>(residuals_).w_ = 1e2; // PositionUpdate
    std::get<7>(residuals_).w_ = 1e4; // GravityCorrectionPrediction
    std::get<8>(residuals_).w_ = 1e2; // AttitudeUpdate
  }
  virtual ~ImuGpsFilter(){}
  virtual void Init(TimePoint t){
    if(GetMinMaxTime() != TimePoint::min()){
      TSIF_LOG("Initializing state at t = " << Print(t));
      startTime_ = t;
      time_ = t;
      state_.SetIdentity();
      state_.Get<1>() = Quat(0.7,0,0,sqrt(1-0.7*0.7));
      I_.setIdentity();
      I_.block<3,3>(State::Start(0),State::Start(0)) = Mat3::Identity()*1e2;
      I_.block<3,3>(State::Start(1),State::Start(1)) = Mat3::Identity()*1e0;
      I_.block<3,3>(State::Start(2),State::Start(2)) = Mat3::Identity()*1e2;
      I_.block<3,3>(State::Start(3),State::Start(3)) = Mat3::Identity()*1e2;
      I_.block<3,3>(State::Start(4),State::Start(4)) = Mat3::Identity()*1e0;
      I_.block<3,3>(State::Start(5),State::Start(5)) = Mat3::Identity()*1e0;
      is_initialized_ = true;
    }
  }
  virtual void ComputeLinearizationPoint(const TimePoint& t){
    curLinState_ = state_;
    double dt = toSec(t-time_);
    std::shared_ptr<const MeasAcc> acc = std::get<2>(timelines_).Get(t);
    std::shared_ptr<const MeasGyr> gyr = std::get<3>(timelines_).Get(t);
    curLinState_.Get<0>() = state_.Get<0>() + dt*state_.Get<1>().toRotationMatrix()*state_.Get<2>();
    curLinState_.Get<1>() = state_.Get<1>()*Exp(dt*state_.Get<3>());
    curLinState_.Get<2>() = (Mat3::Identity() - SSM(dt*state_.template Get<3>()))*state_.template Get<2>()
            + dt*(acc->GetAcc()-state_.template Get<4>()+state_.template Get<1>().inverse().toRotationMatrix()*std::get<2>(residuals_).g_);
    curLinState_.Get<3>() = gyr->GetGyr()-state_.template Get<5>();
  }
};

} // namespace tsif

#endif  // TSIF_IMU_GPS_H_
