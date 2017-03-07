#ifndef TSIF_VIO_H_
#define TSIF_VIO_H_

#include <opencv2/opencv.hpp>

#include "tsif/residuals/accelerometer_prediction.h"
#include "tsif/residuals/attitude_findif.h"
#include "tsif/residuals/bearing_findif.h"
#include "tsif/residuals/gyroscope_update.h"
#include "tsif/residuals/image_update.h"
#include "tsif/residuals/position_findif.h"
#include "tsif/residuals/random_walk.h"
#include "tsif/common.h"
#include "tsif/filter.h"

namespace tsif{

template<int N>
class MeasImg: public ElementVector<>{
 public:
  MeasImg(){}
  MeasImg(const std::string& filename){
    img_ = cv::imread(filename, 1);
  }
  cv::Mat img_;
  mutable std::array<UnitVector,N> n_;
  const UnitVector& GetBea(int i) const{
    return n_[i];
  }
  void SetBea(int i, const UnitVector& n) const{
    n_[i] = n;
  }
  std::vector<cv::KeyPoint> keyPoints_;
};

class VioFilter: public Filter<PositionFindif<0,0,2,1>,
                               AttitudeFindif<0,1,3>,
                               AccelerometerPrediction<0,2,1,3,4>,
                               GyroscopeUpdate<0,3,5>,
                               RandomWalk<Element<Vec3,4>>,
                               RandomWalk<Element<Vec3,5>>,
                               ImageUpdate<0,6,10,MeasImg<10>>,
                               RandomWalk<Element<std::array<UnitVector,10>,6>>,
                               RandomWalk<Element<std::array<double,10>,7>>>{
 public:
  VioFilter(){
    max_iter_ = 10;
    th_iter_ = 0.1;
    include_max_ = true;
    std::get<0>(residuals_).w_ = 1e4; // PositionFindif
    std::get<1>(residuals_).w_ = 1e4; // AttitudeFindif
    std::get<2>(residuals_).w_ = 1e2; // AccelerometerPrediction
    std::get<3>(residuals_).w_ = 1/(2e-3); // GyroscopeUpdate
    std::get<4>(residuals_).w_ = 1e4; // AccelerometerBiasPrediction
    std::get<5>(residuals_).w_ = 1e4; // GyroscopeBiasPrediction
    for(auto i : desc_){
      i = -1;
    }

    BearingFindif<0,0,1,2,3,2> test;
    test.JacPreTest(1e-6,1e-8);
    test.JacCurTest(1e-6,1e-8);
  }
  virtual ~VioFilter(){}
  virtual void Init(TimePoint t){
    if(GetMinMaxTime() != TimePoint::min()){
      TSIF_LOG("Initializing state at t = " << Print(t));
      startTime_ = t;
      time_ = t;
      state_.SetIdentity();
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
  void InitializeNew(std::shared_ptr<const MeasImg<10>> m){
    cv::ORB orb(20);
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat desc;
    orb.detect(m->img_, keyPoints);
    for(auto kp : keyPoints){
      std::cout << kp.pt << "\t" << kp.response << "\t" << kp.octave << std::endl;
    }
    orb.compute(m->img_, keyPoints, desc);
    cv::Mat drawImg;
    cv::drawKeypoints(m->img_,keyPoints,drawImg);

    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Display Image", drawImg);
    cv::waitKey(0);

    // Detect features (mask current)
    // Bucket
    //  - create buckets (add up scores)
    //  - take best from each bucket, starting with highest score bucket
    //  - for each initalize new bearing vector and distance
  }
  void RemoveBad(){

  }
  virtual void PreProcess(){
    for(int i=0;i<10;i++){
      std::get<6>(residuals_).active_[i] = false;
    }
    for(auto kp : std::get<6>(residuals_).meas_->keyPoints_){
      for(int i=0;i<10;i++){
        if(kp.class_id == desc_[i]){
          Vec3 vec;
          cam_.PixelToBearing(Vec<2>(kp.pt.x,kp.pt.y),vec);
          std::get<6>(residuals_).meas_->SetBea(i,UnitVector(vec));
          std::get<6>(residuals_).active_[i] = true;
        }
      }
    }
  };
  virtual void PostProcess(){
    // Remove if not tracked for more then 10 frames
    for(int i=0;i<10;i++){
      if(std::get<6>(residuals_).active_[i]){
        count_[i]++;
      } else {
        count_[i] = 0;
      }
      if(count_[i] > 10){
        desc_[i] = -1;
      }
    }

    // Add new landmarks
    for(int i=0;i<10;i++){
      if(desc_[i] == -1){
        ;
      }
    }
  };
  void AddNewLandmark(int i,const UnitVector& n, int desc){
    count_[i] = 0;
    desc_[i] = desc;
    state_.Get<6>()[i] = n;
    state_.Get<7>()[i] = 1;
    // TODO: covariance
  }
 private:
   std::array<int,10> desc_;
   std::array<int,10> count_;
   Camera cam_;
  // Status [N]
};

} // namespace tsif

#endif  // TSIF_VIO_H_
