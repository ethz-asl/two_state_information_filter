#ifndef TSIF_VIO_H_
#define TSIF_VIO_H_

#include <opencv2/opencv.hpp>

#include "tsif/residuals/accelerometer_prediction.h"
#include "tsif/residuals/attitude_findif.h"
#include "tsif/residuals/bearing_findif.h"
#include "tsif/residuals/distance_findif.h"
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
  MeasImg(){
    for(int i=0;i<N;i++){
      n_.at(i).SetRandom();
    }
  }
  MeasImg(const std::string& filename){
    img_ = cv::imread(filename, 1);
    for(int i=0;i<N;i++){
      n_.at(i).SetRandom();
    }
  }
  cv::Mat img_;
  mutable std::array<UnitVector,N> n_;
  const UnitVector& GetBea(int i) const{
    return n_.at(i);
  }
  void SetBea(int i, const UnitVector& n) const{
    n_.at(i) = n;
  }
  std::vector<cv::KeyPoint> keyPoints_;
};

template<int N>
using VioFilterBase = Filter<PositionFindif<0,0,2,1>,
                             AttitudeFindif<0,1,3>,
                             AccelerometerPrediction<0,2,1,3,4>,
                             GyroscopeUpdate<0,3,5>,
                             RandomWalk<Element<Vec3,4>>,
                             RandomWalk<Element<Vec3,5>>,
                             ImageUpdate<0,6,N,MeasImg<N>>,
                             BearingFindif<0,6,7,2,3,N>,
                             DistanceFindif<0,6,7,2,N>>;

template<int N>
class VioFilter: public VioFilterBase<N> {
  typedef VioFilterBase<N> Base;
  using Base::max_iter_;
  using Base::th_iter_;
  using Base::include_max_;
  using Base::residuals_;
  using Base::startTime_;
  using Base::time_;
  using Base::state_;
  using Base::I_;
  using Base::is_initialized_;
  using Base::curLinState_;
  using Base::timelines_;
  using Base::GetMinMaxTime;
  using typename Base::State;
 public:
  VioFilter(){
    max_iter_ = 10;
    th_iter_ = 0.1;
    include_max_ = false;
    std::get<0>(residuals_).w_ = 1e4; // PositionFindif
    std::get<1>(residuals_).w_ = 1e4; // AttitudeFindif
    std::get<2>(residuals_).w_ = 1e2; // AccelerometerPrediction
    std::get<3>(residuals_).w_ = 1/(2e-3); // GyroscopeUpdate
    std::get<4>(residuals_).w_ = 1e4; // AccelerometerBiasPrediction
    std::get<5>(residuals_).w_ = 1e4; // GyroscopeBiasPrediction
    std::get<6>(residuals_).w_ = 1e2; // Image Update
    std::get<7>(residuals_).w_ = 1e2; // BearingPrediction
    std::get<8>(residuals_).w_ = 1e2; // DistancePrediction
    for(auto& i : desc_){
      i = -1;
    }
  }
  virtual ~VioFilter(){}
  virtual void Init(TimePoint t){
    if(GetMinMaxTime() != TimePoint::min()){
      TSIF_LOG("Initializing state at t = " << Print(t));
      startTime_ = t;
      time_ = t;
      state_.SetIdentity();
      I_.setIdentity();
      I_.template block<3,3>(State::Start(0),State::Start(0)) = Mat3::Identity()*1e2;
      I_.template block<3,3>(State::Start(1),State::Start(1)) = Mat3::Identity()*1e0;
      I_.template block<3,3>(State::Start(2),State::Start(2)) = Mat3::Identity()*1e2;
      I_.template block<3,3>(State::Start(3),State::Start(3)) = Mat3::Identity()*1e2;
      I_.template block<3,3>(State::Start(4),State::Start(4)) = Mat3::Identity()*1e0;
      I_.template block<3,3>(State::Start(5),State::Start(5)) = Mat3::Identity()*1e0;
      is_initialized_ = true;
    }
  }
  virtual void ComputeLinearizationPoint(const TimePoint& t){
    curLinState_ = state_;
    double dt = toSec(t-time_);
    std::shared_ptr<const MeasAcc> acc = std::get<2>(timelines_).Get(t);
    std::shared_ptr<const MeasGyr> gyr = std::get<3>(timelines_).Get(t);
    curLinState_.template Get<0>() = state_.template Get<0>() + dt*state_.template Get<1>().toRotationMatrix()*state_.template Get<2>();
    curLinState_.template Get<1>() = state_.template Get<1>()*Exp(dt*state_.template Get<3>());
    curLinState_.template Get<2>() = (Mat3::Identity() - SSM(dt*state_.template Get<3>()))*state_.template Get<2>()
            + dt*(acc->GetAcc()-state_.template Get<4>()+state_.template Get<1>().inverse().toRotationMatrix()*std::get<2>(residuals_).g_);
    curLinState_.template Get<3>() = gyr->GetGyr()-state_.template Get<5>();
  }
  void InitializeNew(std::shared_ptr<const MeasImg<N>> m){
    cv::ORB orb(20);
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat desc;
    orb.detect(m->img_, keyPoints);
    for(const auto& kp : keyPoints){
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
    drawImg_ = cv::Mat::zeros(480, 752, CV_8U);
    for(int i=0;i<N;i++){
      std::get<6>(residuals_).active_[i] = false;
    }
    cv::drawKeypoints(drawImg_,std::get<6>(residuals_).meas_->keyPoints_,drawImg_,cv::Scalar(255,50,50));
    std::vector<cv::KeyPoint> matches;
    for(const auto& kp : std::get<6>(residuals_).meas_->keyPoints_){
      for(int i=0;i<N;i++){
        if(kp.class_id == desc_[i]){
          TSIF_LOGW("Found measurement for landmark " << i);
          Vec<2> pt_pred;
          cam_.BearingToPixel(state_.template Get<6>()[i].GetVec(),pt_pred);
          matches.emplace_back(pt_pred(0),pt_pred(1),0,0,0,0,desc_[i]);
          Vec3 vec;
          cam_.PixelToBearing(Vec<2>(kp.pt.x,kp.pt.y),vec);
          std::get<6>(residuals_).meas_->SetBea(i,UnitVector(vec));
          std::get<6>(residuals_).active_[i] = true;
        }
      }
    }
    cv::drawKeypoints(drawImg_,matches,drawImg_,cv::Scalar(50,50,255));
  };
  virtual void PostProcess(){
    // Remove if not tracked for more then 10 frames
    for(int i=0;i<N;i++){
      if(std::get<6>(residuals_).active_[i] || desc_[i] == -1){
        count_[i] = 0;
      } else {
        count_[i]++;
      }
      if(count_[i] > 10){
        TSIF_LOGW("Removing landmark " << i);
        desc_[i] = -1;
      }
    }

    // Add new landmarks
    for(int i=0;i<N;i++){
      if(desc_[i] == -1){
        for(const auto& kp : std::get<6>(residuals_).meas_->keyPoints_){
          bool isAlreadyUsed = false;
          for(int j=0;j<N;j++){
            if(kp.class_id == desc_[j]){
              isAlreadyUsed = true;
              break;
            }
          }
          if(!isAlreadyUsed){
            TSIF_LOGW("Adding landmark " << i << " with desc " << kp.class_id);
            Vec3 vec;
            cam_.PixelToBearing(Vec<2>(kp.pt.x,kp.pt.y),vec);
            AddNewLandmark(i,UnitVector(vec),kp.class_id,1);
            break;
          }
        }
      }
    }
    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Display Image", drawImg_);
    cv::waitKey(1);
  };
  void AddNewLandmark(int i,const UnitVector& n, int desc, double indDis){
    count_[i] = 0;
    desc_[i] = desc;
    state_.template Get<6>()[i] = n;
    state_.template Get<7>()[i] = indDis;
    I_.template block<2,State::Dim()>(State::Start(6)+2*i,0).setZero();
    I_.template block<1,State::Dim()>(State::Start(7)+1*i,0).setZero();
    I_.template block<State::Dim(),2>(0,State::Start(6)+2*i).setZero();
    I_.template block<State::Dim(),1>(0,State::Start(7)+1*i).setZero();
    I_.template block<2,2>(State::Start(6)+2*i,State::Start(6)+2*i) = Mat<2>::Identity()*1e2;
    I_.template block<1,1>(State::Start(7)+1*i,State::Start(7)+1*i) = Mat<1>::Identity()*1e-4;
  }
 private:
   std::array<int,N> desc_;
   std::array<int,N> count_; // Status
   Camera cam_;
   cv::Mat drawImg_;
};

} // namespace tsif

#endif  // TSIF_VIO_H_
