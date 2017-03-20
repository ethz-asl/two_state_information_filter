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
    isSim_ = false;
  }
  MeasImg(const std::string& filename){
    img_ = cv::imread(filename, 0);
    for(int i=0;i<N;i++){
      n_.at(i).SetRandom();
    }
    isSim_ = false;
  }
  const UnitVector& GetBea(int i) const{
    return n_.at(i);
  }
  void SetBea(int i, const UnitVector& n) const{
    n_.at(i) = n;
  }
  cv::Mat img_;
  mutable std::array<UnitVector,N> n_;
  mutable std::vector<cv::KeyPoint> keyPoints_;
  bool isSim_;
};

class LandmarkData{
 public:
  LandmarkData(){
    id_ = -1;
  }
  Vec<2> prePoint_;
  Mat<2> preCov_;
  cv::Mat desc_;
  double trackingQuality_;
  std::vector<std::vector<cv::DMatch>> matches_;
  int count_;
  double sigmaAngle_;
  double sigma0_;
  double sigma1_;
  int id_;

  void SetPrediction(const Camera& cam, const UnitVector& vec, MatCRef<2> P){
    UnitVector vec0, vec1;
    Vec<2> pt_pred0, pt_pred1;
    cam.BearingToPixel(vec.GetVec(),prePoint_);
    Eigen::EigenSolver<Mat<2>> ES(P);
    vec.Boxplus(ES.eigenvectors().col(0).real()*sqrt(ES.eigenvalues()(0).real()),vec0);
    vec.Boxplus(ES.eigenvectors().col(1).real()*sqrt(ES.eigenvalues()(1).real()),vec1);
    cam.BearingToPixel(vec0.GetVec(),pt_pred0);
    cam.BearingToPixel(vec1.GetVec(),pt_pred1);
    Mat<2> U;
    U.col(0) = pt_pred0-prePoint_;
    U.col(1) = pt_pred1-prePoint_;
    preCov_ = (U*U.transpose());
    ES.compute(preCov_);
    sigmaAngle_ = std::atan2(ES.eigenvectors()(1,0).real(),ES.eigenvectors()(0,0).real());
    sigma0_ = sqrt(ES.eigenvalues()(0).real());
    sigma1_ = sqrt(ES.eigenvalues()(1).real());
  }
};

template<int N>
using VioFilterBase = Filter<PositionFindif<0,0,2,1>,
                             AttitudeFindif<0,1,3>,
                             AccelerometerPrediction<0,2,1,3,4>,
                             GyroscopeUpdate<0,3,5>,
                             RandomWalk<Element<Vec3,4>>,
                             RandomWalk<Element<Vec3,5>>,
                             ImageUpdate<0,6,N,MeasImg<N>>,
                             BearingFindif<0,6,7,2,3,8,9,N>,
                             DistanceFindif<0,6,7,2,3,8,9,N>,
                             RandomWalk<Element<Vec3,8>>,
                             RandomWalk<Element<Quat,9>>>;

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
  VioFilter(const std::string& optionFile){
    optionFile_ = optionFile;
    max_iter_ = tsif::OptionLoader::Instance().Get<int>(optionFile,"max_iter");
    th_iter_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"th_iter");
    include_max_ = false;
    std::get<0>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_posfd");
    std::get<1>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_attfd");
    std::get<2>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_accpre");
    std::get<3>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_gyrupd");
    std::get<4>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_accbias");
    std::get<5>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_gyrbias");
    std::get<6>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_imgupd");
    std::get<7>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_beapre");
    std::get<8>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_dispre");
    std::get<9>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_veppre");
    std::get<10>(residuals_).w_ = tsif::OptionLoader::Instance().Get<double>(optionFile,"w_veapre");
    drawAdding_ = tsif::OptionLoader::Instance().Get<int>(optionFile_,"draw_adding");
    drawTracking_ = tsif::OptionLoader::Instance().Get<int>(optionFile_,"draw_tracking");
    doDraw_ = tsif::OptionLoader::Instance().Get<int>(optionFile_,"do_draw");
    singleMatchOnly_ = tsif::OptionLoader::Instance().Get<int>(optionFile_,"single_match_only");
  }
  virtual ~VioFilter(){}
  virtual void Init(TimePoint t){
    if(GetMinMaxTime() != TimePoint::min()){
      TSIF_LOG("Initializing state at t = " << Print(t));
      startTime_ = t;
      time_ = t;
      state_.SetIdentity();
      state_.template Get<1>() = Quat::FromTwoVectors(std::get<2>(timelines_).Get(std::get<2>(timelines_).GetFirstTime())->GetAcc(),Vec3(0,0,1));
      state_.template Get<8>() = tsif::OptionLoader::Instance().Get<Vec3>(optionFile_,"init_vep");
      state_.template Get<9>() = tsif::OptionLoader::Instance().Get<Quat>(optionFile_,"init_vea");
      I_.setIdentity();
      I_.template block<3,3>(State::Start(0),State::Start(0)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_pos"),2);
      I_.template block<3,3>(State::Start(1),State::Start(1)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_att"),2);
      I_.template block<3,3>(State::Start(2),State::Start(2)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_vel"),2);
      I_.template block<3,3>(State::Start(3),State::Start(3)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_ror"),2);
      I_.template block<3,3>(State::Start(4),State::Start(4)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_acb"),2);
      I_.template block<3,3>(State::Start(5),State::Start(5)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_gyb"),2);
      I_.template block<3,3>(State::Start(8),State::Start(8)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_vep"),2);
      I_.template block<3,3>(State::Start(9),State::Start(9)) = Mat3::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_vea"),2);
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

    const Mat3 C_VI = state_.template Get<9>().toRotationMatrix();
    const Vec3 ror = C_VI*state_.template Get<3>();
    const Vec3 vel = C_VI*(state_.template Get<2>() + state_.template Get<3>().cross(state_.template Get<8>()));
    for(int i=0;i<N;i++){
      const UnitVector& bea = state_.template Get<6>()[i];
      const Vec3 beaVec = bea.GetVec();
      const Mat<3,2> beaN = bea.GetN();
      const double& invDis = state_.template Get<7>()[i];
      const Vec<2> dn = -beaN.transpose()*(ror + invDis * beaVec.cross(vel));
      bea.Boxplus(dn*dt,curLinState_.template Get<6>()[i]);
      curLinState_.template Get<7>()[i] = invDis + dt*beaVec.dot(vel)*invDis*invDis;
    }
  }
  virtual void PreProcess(){
    const auto& m = std::get<6>(residuals_).meas_;
    for(int i=0;i<N;i++){
      std::get<6>(residuals_).active_[i] = false;
    }

    Timer timer;
    if(doDraw_){
      if(!m->isSim_){
        cv::cvtColor(m->img_,drawImg_,CV_GRAY2BGR);
      } else {
        drawImg_ = cv::Mat::zeros(480, 752, CV_8U);
        cv::drawKeypoints(drawImg_,m->keyPoints_,drawImg_,cv::Scalar(255,50,50));
      }
    }
    std::cout << "Img: " << 1000*timer.GetIncr() << std::endl;

    // Compute covariance of landmark predition // TODO: make more efficient
    MatX P = I_.inverse();
    MatX JBearingPrediction(2*N,State::Dim());
    JBearingPrediction.setZero();
    std::get<7>(residuals_).JacPreCustom(JBearingPrediction,state_,curLinState_,true);
    MatX Pbearing = (JBearingPrediction*P*JBearingPrediction.transpose()) + MatX(Vec<2*N>::Ones().asDiagonal()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"w_imgupd"),-2));
    std::cout << "Cov: " << 1000*timer.GetIncr() << std::endl;

    // Detect and compute
    if(!m->isSim_){
      mask_ = cv::Mat::zeros(m->img_.size(), CV_8U);
      std::vector<cv::KeyPoint> keyPoints;
      cv::Mat desc;
      for(int i=0;i<N;i++){
        if(!l_.at(i).desc_.empty()){
          // Extract pixel coordinates and corresponding uncertainty
          l_.at(i).SetPrediction(cam_,curLinState_.template Get<6>()[i],Pbearing.block<2,2>(2*i,2*i));
          cv::ellipse(mask_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Size(l_.at(i).sigma0_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance"),l_.at(i).sigma1_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance")), l_.at(i).sigmaAngle_/M_PI*180, 0, 360, cv::Scalar(255), -1);
        }
      }
      std::cout << "Mask: " << 1000*timer.GetIncr() << std::endl;
      cv::ORB orb(tsif::OptionLoader::Instance().Get<int>(optionFile_,"num_candidates_matching")*NumLandmarks(),2,1);
      orb.detect(m->img_,keyPoints,mask_);
      std::cout << "Detection: " << 1000*timer.GetIncr() << std::endl;
      orb.compute(m->img_,keyPoints,desc);
      std::cout << "Compute: " << 1000*timer.GetIncr() << std::endl;
      if(doDraw_ && drawTracking_){
        cv::drawKeypoints(drawImg_,keyPoints,drawImg_,cv::Scalar(255,0,0));
      }

      for(int i=0;i<N;i++){
        if(!l_.at(i).desc_.empty()){
          // Match landmarks
          cv::BFMatcher matcher(cv::NORM_HAMMING);
          l_.at(i).matches_.clear();
          matcher.radiusMatch(l_.at(i).desc_,desc,l_.at(i).matches_,tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_desc_distance"));
          int best_match = -1;
          double best_score = 0;
          int count = 0;
          for(int j=0;j<l_.at(i).matches_[0].size();j++){
            const Vec<2> pt_meas(keyPoints[l_.at(i).matches_[0][j].trainIdx].pt.x,keyPoints[l_.at(i).matches_[0][j].trainIdx].pt.y);
            const double geomScore = std::sqrt(((l_.at(i).prePoint_-pt_meas).transpose()*l_.at(i).preCov_.inverse()*(l_.at(i).prePoint_-pt_meas))(0))/tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance");
            const double descScore = l_.at(i).matches_[0][j].distance/tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_desc_distance");
            if(geomScore + descScore < 1){
              count ++;
              if(best_match == -1 || geomScore + descScore < best_score){
                best_match = j;
                best_score = geomScore + descScore;
              }
            }
          }
          if(best_match != -1 && (!singleMatchOnly_ || count == 1)){
            TSIF_LOGW("Found measurement for landmark " << i);
            Vec<2> pt_meas(keyPoints[l_.at(i).matches_[0][best_match].trainIdx].pt.x,keyPoints[l_.at(i).matches_[0][best_match].trainIdx].pt.y);
            Vec3 vec;
            cam_.PixelToBearing(pt_meas,vec);
            m->SetBea(i,UnitVector(vec));
            std::get<6>(residuals_).active_[i] = true;
            if(doDraw_ && drawTracking_){
              cv::ellipse(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Size(l_.at(i).sigma0_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance"),l_.at(i).sigma1_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance")), l_.at(i).sigmaAngle_/M_PI*180, 0, 360, cv::Scalar(0,255,0));
              cv::line(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Point2f(pt_meas(0),pt_meas(1)), cv::Scalar(0,255,0));
            }
          } else {
            if(doDraw_ && drawTracking_){
              if(count == 0){
                cv::ellipse(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Size(l_.at(i).sigma0_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance"),l_.at(i).sigma1_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance")), l_.at(i).sigmaAngle_/M_PI*180, 0, 360, cv::Scalar(0,0,255));
              } else {
                cv::ellipse(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Size(l_.at(i).sigma0_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance"),l_.at(i).sigma1_*tsif::OptionLoader::Instance().Get<float>(optionFile_,"match_geom_distance")), l_.at(i).sigmaAngle_/M_PI*180, 0, 360, cv::Scalar(255,0,255));
              }
            }
          }
        }
      }
    } else {
      for(int i=0;i<N;i++){
        std::vector<cv::KeyPoint> matches;
        for(const auto& kp : std::get<6>(residuals_).meas_->keyPoints_){
          for(int i=0;i<N;i++){
            if(kp.class_id == l_.at(i).id_){
              TSIF_LOGW("Found measurement for landmark " << i);
              l_.at(i).SetPrediction(cam_,curLinState_.template Get<6>()[i],Pbearing.block<2,2>(2*i,2*i)); // TODO: change to curLinState_
              matches.emplace_back(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1),0,0,0,0,l_.at(i).id_);
              Vec3 vec;
              cam_.PixelToBearing(Vec<2>(kp.pt.x,kp.pt.y),vec);
              m->SetBea(i,UnitVector(vec));
              std::get<6>(residuals_).active_[i] = true;
              if(doDraw_ && drawTracking_){
                cv::ellipse(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), cv::Size(l_.at(i).sigma0_,l_.at(i).sigma1_), l_.at(i).sigmaAngle_/M_PI*180, 0, 360, cv::Scalar(0,255,0));
                cv::line(drawImg_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), kp.pt, cv::Scalar(0,255,0));
              }
            }
          }
        }
      }
    }
    std::cout << "Matching: " << 1000*timer.GetIncr() << std::endl;

    if(doDraw_){
      DrawVirtualHorizon(drawImg_,state_.template Get<9>()*state_.template Get<1>().inverse());
    }
    std::cout << "Draw Horizon: " << 1000*timer.GetIncr() << std::endl;
    std::cout << "=== Preprocess: " << 1000*timer.GetFull() << std::endl;
  };
  void DrawVirtualHorizon(cv::Mat& img,Quat att){
    cv::rectangle(img,cv::Point2f(0,0),cv::Point2f(82,92),cv::Scalar(50,50,50),-1,8,0);
    cv::rectangle(img,cv::Point2f(0,0),cv::Point2f(80,90),cv::Scalar(100,100,100),-1,8,0);
    cv::Point2f rollCenter = cv::Point2f(40,40);
    cv::Scalar rollColor1(50,50,50);
    cv::Scalar rollColor2(200,200,200);
    cv::Scalar rollColor3(120,120,120);
    cv::circle(img,rollCenter,32,rollColor1,-1,8,0);
    cv::circle(img,rollCenter,30,rollColor2,-1,8,0);
    Eigen::Vector3d Vg = att.toRotationMatrix()*Eigen::Vector3d(0,0,-1);
    double roll = atan2(Vg(1),Vg(0))-0.5*M_PI;
    double pitch = acos(Vg.dot(Eigen::Vector3d(0,0,1)))-0.5*M_PI;
    double pixelFor10Pitch = 5.0;
    double pitchOffsetAngle = -asin(pitch/M_PI*180.0/10.0*pixelFor10Pitch/30.0);
    cv::Point2f rollVector1 = 30*cv::Point2f(cos(roll),sin(roll));
    cv::Point2f rollVector2 = cv::Point2f(25,0);
    cv::Point2f rollVector3 = cv::Point2f(10,0);
    std::vector<cv::Point> pts;
    cv::ellipse2Poly(rollCenter,cv::Size(30,30),0,(roll-pitchOffsetAngle)/M_PI*180,(roll+pitchOffsetAngle)/M_PI*180+180,1,pts);
    cv::Point *points;
    points = &pts[0];
    int nbtab = pts.size();
    cv::fillPoly(img,(const cv::Point**)&points,&nbtab,1,rollColor3);
    cv::line(img,rollCenter+rollVector2,rollCenter+rollVector3,rollColor1, 2);
    cv::line(img,rollCenter-rollVector2,rollCenter-rollVector3,rollColor1, 2);
    cv::ellipse(img,rollCenter,cv::Size(10,10),0,0,180,rollColor1,2,8,0);
    cv::circle(img,rollCenter,2,rollColor1,-1,8,0);
  }
  virtual void PostProcess(){
    Timer timer;
    const auto& m = std::get<6>(residuals_).meas_;
    // Remove if not tracked for more then prune_count frames
    for(int i=0;i<N;i++){
      if(std::get<6>(residuals_).active_[i] || (!m->isSim_ && l_.at(i).desc_.empty()) || (m->isSim_ && l_.at(i).id_ == -1)){
        l_.at(i).count_ = 0;
      } else {
        l_.at(i).count_++;
      }
      if(l_.at(i).count_ > tsif::OptionLoader::Instance().Get<int>(optionFile_,"prune_count") || state_.template Get<7>()[i] <= 0){
        TSIF_LOGW("Removing landmark " << i);
        l_.at(i).desc_ = cv::Mat();
        l_.at(i).id_ = -1;
      }
    }
    std::cout << "Removing: " << 1000*timer.GetIncr() << std::endl;

    // Add new landmarks
    if(N - NumLandmarks() > tsif::OptionLoader::Instance().Get<int>(optionFile_,"start_add")){

      // Create mask
      std::vector<cv::KeyPoint> keyPoints;
      cv::Mat desc;
      mask_ = cv::Mat::ones(m->img_.size(), CV_8U); // TODO: adapt size to sim
      for(int i=0;i<N;i++){
        if(!l_.at(i).desc_.empty()){
          cv::circle(mask_, cv::Point2f(l_.at(i).prePoint_(0),l_.at(i).prePoint_(1)), tsif::OptionLoader::Instance().Get<int>(optionFile_,"neighbor_suppression"), cv::Scalar(0), -1); // TODO: param
        }
      }
      std::cout << "Mask: " << 1000*timer.GetIncr() << std::endl;

      if(!m->isSim_){
        cv::ORB orb(tsif::OptionLoader::Instance().Get<int>(optionFile_,"num_candidates_add"),2,1);
        orb.detect(m->img_,keyPoints,mask_);
        std::cout << "Detection: " << 1000*timer.GetIncr() << std::endl;
        orb.compute(m->img_,keyPoints,desc);
        std::cout << "Compute: " << 1000*timer.GetIncr() << std::endl;
      } else {
        keyPoints = m->keyPoints_;
      }

      if(doDraw_ && drawAdding_){
        cv::drawKeypoints(drawImg_,keyPoints,drawImg_,cv::Scalar(255,50,50));
      }
      std::cout << "Drawing Candidates: " << 1000*timer.GetIncr() << std::endl;
      const int B = tsif::OptionLoader::Instance().Get<int>(optionFile_,"bucket_count");
      int bucketsCandidates[B*B];
      for(int i=0;i<B;i++){
        for(int j=0;j<B;j++){
          bucketsCandidates[j*B+i] = -1;
        }
        if(doDraw_ && drawAdding_){
          cv::line(drawImg_, cv::Point(752*i/B,0), cv::Point(752*i/B,480), cv::Scalar(0,0,255));
          cv::line(drawImg_, cv::Point(0,480*i/B), cv::Point(752,480*i/B), cv::Scalar(0,0,255));
        }
      }

      // Compute buckets
      for(int i=0;i<keyPoints.size();i++){
        const auto& kp = keyPoints[i];
        int x_ind = std::floor(kp.pt.x/752*B);
        int y_ind = std::floor(kp.pt.y/480*B);
        if(x_ind >= 0 && x_ind < B && y_ind >= 0 && y_ind < B){
          if(m->isSim_){
            bool found = false;
            for(int j=0;j<N;j++){
              found = found | l_.at(j).id_ == kp.class_id;
            }
            if(found){
              continue;
            }
          }
          if(bucketsCandidates[y_ind*B+x_ind] == -1 || kp.response > keyPoints[bucketsCandidates[y_ind*B+x_ind]].response){
            bucketsCandidates[y_ind*B+x_ind] = i;
          }
        }
      }
      std::cout << "Bucketing: " << 1000*timer.GetIncr() << std::endl;

      std::vector<cv::KeyPoint> newKeyPoints;
      for(int i=0;i<N;i++){
        if((!m->isSim_ && l_.at(i).desc_.empty()) || (m->isSim_ && l_.at(i).id_ == -1)){
          // Look for non-empty bucket with no landmark // TODO: start with best
          const int start = abs((int)(NormalRandomNumberGenerator::Instance().Get()*1e6));
          for(int j=0;j<B*B;j++){
            if(bucketsCandidates[(j+start)%(B*B)] >= 0){
              const int newKeyPointInd = bucketsCandidates[(j+start)%(B*B)];
              newKeyPoints.push_back(keyPoints[newKeyPointInd]);
              bucketsCandidates[(j+start)%(B*B)] = -1;
              if(!m->isSim_){
                desc.row(newKeyPointInd).copyTo(l_.at(i).desc_);
              } else {
                l_.at(i).id_ = keyPoints[newKeyPointInd].class_id;
              }
              TSIF_LOGW("Adding landmark " << i);
              Vec3 vec;
              cam_.PixelToBearing(Vec<2>(keyPoints[newKeyPointInd].pt.x,keyPoints[newKeyPointInd].pt.y),vec);
              AddNewLandmark(i,UnitVector(vec),tsif::OptionLoader::Instance().Get<double>(optionFile_,"init_dis"));
              break;
            }
          }
        }
      }
      std::cout << "Adding: " << 1000*timer.GetIncr() << std::endl;
      if(doDraw_ && drawAdding_){
        cv::drawKeypoints(drawImg_,newKeyPoints,drawImg_,cv::Scalar(255,50,255));
      }
    }

    if(doDraw_){
      cv::namedWindow("VIO", cv::WINDOW_AUTOSIZE);
      cv::imshow("VIO", drawImg_);
      cv::waitKey(2);
    }
    std::cout << "Drawing: " << 1000*timer.GetIncr() << std::endl;
    std::cout << "=== Postprocess: " << 1000*timer.GetFull() << std::endl;
  };
  void AddNewLandmark(int i,const UnitVector& n, double indDis){
    l_.at(i).count_= 0;
    state_.template Get<6>()[i] = n;
    state_.template Get<7>()[i] = indDis;
    I_.template block<2,State::Dim()>(State::Start(6)+2*i,0).setZero();
    I_.template block<1,State::Dim()>(State::Start(7)+1*i,0).setZero();
    I_.template block<State::Dim(),2>(0,State::Start(6)+2*i).setZero();
    I_.template block<State::Dim(),1>(0,State::Start(7)+1*i).setZero();
    I_.template block<2,2>(State::Start(6)+2*i,State::Start(6)+2*i) = Mat<2>::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_bea"),2);
    I_.template block<1,1>(State::Start(7)+1*i,State::Start(7)+1*i) = Mat<1>::Identity()*pow(tsif::OptionLoader::Instance().Get<double>(optionFile_,"initstd_dis"),2);
  }
  int NumLandmarks(){
    int count = 0;
    for(int i=0;i<N;i++){
      count += (!l_.at(i).desc_.empty() | l_.at(i).id_ != -1);
    }
    return count;
  }
 private:
  Camera cam_;
  cv::Mat drawImg_;
  std::string optionFile_;
  std::array<LandmarkData,N> l_;
  bool drawAdding_;
  bool drawTracking_;
  bool doDraw_;
  bool singleMatchOnly_;
  cv::Mat mask_;
};

} // namespace tsif

#endif  // TSIF_VIO_H_
