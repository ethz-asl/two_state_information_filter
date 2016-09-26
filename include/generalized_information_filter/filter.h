#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/measurement.h"

namespace GIF {

/*! \brief Residual Struct
 *         Contains various object and temporaries associated with a specific residual.
 */
struct ResidualStruct {
  ResidualStruct(const std::shared_ptr<BinaryResidualBase>& res,
                 const std::shared_ptr<ElementVectorDefinition>& stateDefinition,
                 const std::shared_ptr<ElementVectorDefinition>& noiseDefinition,
                 const Duration& maxWaitTime,
                 const Duration& minWaitTime,
                 const std::string& name):
                     mt_(!res->isUnary_, maxWaitTime, minWaitTime),
                     preWrap_(res->PreDefinition(), stateDefinition),
                     curWrap_(res->CurDefinition(), stateDefinition),
                     noiWrap_(res->NoiDefinition(), noiseDefinition),
                     inn_(res->InnDefinition()),
                     innRef_(res->InnDefinition()),
                     name_(name){
    stateDefinition->ExtendWithOtherElementVectorDefinition(*res->PreDefinition());
    stateDefinition->ExtendWithOtherElementVectorDefinition(*res->CurDefinition());
    noiseDefinition->ExtendWithOtherElementVectorDefinition(*res->NoiDefinition());
    res_ = res;
    preWrap_.ComputeMap();
    curWrap_.ComputeMap();
    noiWrap_.ComputeMap();
    innRef_.SetIdentity();
    innDim_ = res->InnDefinition()->GetDim();
    jacPre_.resize(innDim_, res->PreDefinition()->GetDim());
    jacCur_.resize(innDim_, res->CurDefinition()->GetDim());
    jacNoi_.resize(innDim_, res->NoiDefinition()->GetDim());
  }
  ~ResidualStruct() {
  }

  std::shared_ptr<BinaryResidualBase> res_;
  MeasurementTimeline mt_;
  ElementVectorWrapper preWrap_;
  ElementVectorWrapper curWrap_;
  ElementVectorWrapper noiWrap_;
  ElementVector inn_;
  ElementVector innRef_;
  MatX jacPre_;
  MatX jacCur_;
  MatX jacNoi_;
  int innDim_;
  std::string name_;
};

/*! \brief Filter
 *         Handles residuals. Builds state defintion. Contains measurements timelines. Implements
 *         timing logic.
 */
class Filter {
 public:
  Filter(): stateDefinition_(new ElementVectorDefinition()),
            state_(stateDefinition_), curLinState_(stateDefinition_),
            noiseDefinition_(new ElementVectorDefinition()),
            noise_(noiseDefinition_){
    is_initialized_ = false;
  }

  virtual ~Filter() {
  }

  virtual void Init(const TimePoint& t, const ElementVectorBase* initState = nullptr) {
    LOG(INFO) << "Initializing state at t = "
        << std::chrono::system_clock::to_time_t(t) << std::endl;
    startTime_ = t;
    time_ = t;
    if(initState != nullptr){
      state_ = *initState;
    } else {
      state_.SetIdentity();
    }
    cov_.setIdentity();
    is_initialized_ = true;
  }

  void Construct(){
    state_.Construct();
    state_.SetIdentity();
    curLinState_.Construct();
    noise_.Construct();
    noise_.SetIdentity();
    cov_.resize(stateDefinition_->GetDim(), stateDefinition_->GetDim());
  }

  int AddResidual(const std::shared_ptr<BinaryResidualBase>& res,
                  const Duration& maxWaitTime,
                  const Duration& minWaitTime,
                  const std::string& name) {
    residuals_.emplace_back(res, stateDefinition_, noiseDefinition_, maxWaitTime, minWaitTime,name);
    Construct();
    return residuals_.size() - 1;
  }

  std::shared_ptr<ElementVectorDefinition> StateDefinition() const {
    return stateDefinition_;
  }

  std::shared_ptr<ElementVectorDefinition> NoiseDefinition() const {
    return noiseDefinition_;
  }

  TimePoint GetCurrentTimeFromMeasurements() const {
    TimePoint currentTime = TimePoint::min();
    for (int i = 0; i < residuals_.size(); i++) {
      currentTime = std::max(currentTime, residuals_.at(i).mt_.GetLastTime());
    }
    return currentTime;
  }

  TimePoint GetMaxUpdateTime(const TimePoint& currentTime) const {
    TimePoint maxUpdateTime = currentTime;
    for (int i = 0; i < residuals_.size(); i++) {
      maxUpdateTime = std::min(maxUpdateTime,
                               residuals_.at(i).mt_.GetMaximalUpdateTime(currentTime));
    }
    return maxUpdateTime;
  }

  TimePoint GetMinMeasTime() const {
    TimePoint minMeasTime = TimePoint::max();
    for (int i = 0; i < residuals_.size(); i++) {
      minMeasTime = std::min(minMeasTime, residuals_.at(i).mt_.GetFirstTime());
    }
    return minMeasTime;
  }

  void GetMeasurementTimeList(std::set<TimePoint>& times,
                              const TimePoint& maxUpdateTime,
                              const bool includeMax) const {
    for (int i = 0; i < residuals_.size(); i++) {
      if (!residuals_.at(i).res_->isMergeable_) {
        // Add all non-mergeable measurement times
        residuals_.at(i).mt_.GetAllInRange(times, time_, maxUpdateTime);
      } else if (!residuals_.at(i).res_->isSplitable_ && residuals_.at(i).res_->isUnary_) {
        // For the special case of unary and mergeable residuals add the last measurement time
        residuals_.at(i).mt_.GetLastInRange(times, time_, maxUpdateTime);
      }
    }
    if (includeMax && maxUpdateTime > time_) {
      times.insert(maxUpdateTime);
    }
  }

  void SplitAndMergeMeasurements(const std::set<TimePoint>& times) {
    for (int i = 0; i < residuals_.size(); i++) {
      if (residuals_.at(i).res_->isSplitable_ && !residuals_.at(i).res_->isUnary_) {
        // First insert all (splitable + !unary)
        residuals_.at(i).mt_.Split(times, residuals_.at(i).res_.get());
      }
      if (residuals_.at(i).res_->isMergeable_) {
        // Then remove all undesired (mergeable)
        residuals_.at(i).mt_.MergeUndesired(times, residuals_.at(i).res_.get());
      }
    }
  }

  void AddMeasurement(const int i, const std::shared_ptr<const ElementVectorBase>& meas,
                      const TimePoint& t) {
    if(t <= time_){
      LOG(WARNING) << "Adding measurements before current time (will be discarded)" << std::endl;
      return;
    }
    LOG(INFO) << "Adding measurement with ID " << i << " at  t = "
        << std::chrono::system_clock::to_time_t(t) << std::endl;
    if(residuals_.at(i).res_->CheckMeasType(meas)){
      residuals_.at(i).mt_.AddMeasurement(meas, t);
    } else {
      LOG(ERROR) << "Passing wrong measurement type";
    }
  }

  void PrintMeasurementTimelines(const TimePoint& start, int startOffset, double resolution) {
    std::ostringstream out;
    for (int i = 0; i < startOffset; i++) {
      out << " ";
    }
    out << "|";
    LOG(INFO) << out;
    for (int i = 0; i < residuals_.size(); i++) {
      LOG(INFO) << residuals_.at(i).mt_.Print(start, startOffset, resolution);
    }
  }

  void AttemptInitialization(){
    bool check = true;
    for (int i = 0; i < residuals_.size(); i++) {
      if(residuals_.at(i).mt_.GetFirstTime() == TimePoint::max()
         && !residuals_.at(i).res_->isUnary_){
        check = false;
      }
    }
    if(check){
      Init(GetMinMeasTime());
    }
  }

  void Update() {
    // Initialize if possible
    if(!is_initialized_ && GetMinMeasTime() != TimePoint::max()){
      AttemptInitialization();
    }

    if(is_initialized_){
      for (int i = 0; i < residuals_.size(); i++) {
        // Remove outdated
        residuals_.at(i).mt_.RemoveOutdated(time_);
        // Check time
        if(residuals_.at(i).mt_.GetLastProcessedTime() > time_){
          LOG(ERROR) << "Last processed time is in future, this should not happen!" << std::endl;
        }
      }
      PrintMeasurementTimelines(time_, 20, 0.001);
      LOG(INFO) << "stateTime:\t" << toSec(time_ - startTime_);
      TimePoint currentTime = GetCurrentTimeFromMeasurements();
      LOG(INFO) << "currentTime:\t" << toSec(currentTime - startTime_);
      TimePoint maxUpdateTime = GetMaxUpdateTime(currentTime);
      LOG(INFO) << "maxUpdateTime:\t" << toSec(maxUpdateTime - startTime_);
      std::set<TimePoint> times;
      GetMeasurementTimeList(times, maxUpdateTime, false);
      std::ostringstream out;
      out << "updateTimes:\t";
      for (const auto& t : times) {
        out << toSec(t - startTime_) << "\t";
      }
      LOG(INFO) << out;
      SplitAndMergeMeasurements(times);
      PrintMeasurementTimelines(time_, 20, 0.001);

      // Carry out updates
      for (const auto& t : times) {
        MakeUpdateStep(t);
      }
    }
  }

  void MakeUpdateStep(const TimePoint& t) {
    // Compute linearisation point
    curLinState_ = state_;

    // Check available measurements and prepare residuals
    int innDim = 0;
    std::vector<bool> hasMeas(residuals_.size(), false);
    std::shared_ptr<const ElementVectorBase> meas;
    for (int i = 0; i < residuals_.size(); i++) {
      hasMeas.at(i) = residuals_.at(i).mt_.GetMeasurement(t, meas);
      if (hasMeas.at(i)) {
        residuals_.at(i).res_->SetDt(toSec(t-time_));
        residuals_.at(i).res_->SetMeas(meas);
        PreProcess(i);
        innDim += residuals_.at(i).res_->InnDefinition()->GetDim();
      }
    }

    // Temporaries
    VecX y(innDim);
    MatX JacPre(innDim, stateDefinition_->GetDim());
    MatX JacCur(innDim, stateDefinition_->GetDim());
    MatX JacNoi(innDim, noiseDefinition_->GetDim());
    JacPre.setZero();
    JacCur.setZero();
    JacNoi.setZero();
    MatX R(noiseDefinition_->GetDim(), noiseDefinition_->GetDim());
    MatX Winv(innDim, innDim);
    R.setIdentity();
    Winv.setZero();

    // Evaluate residuals and Jacobians
    int count = 0;
    for (int i = 0; i < residuals_.size(); i++) {
      if (hasMeas.at(i)) {
        residuals_.at(i).preWrap_.SetElementVector(&state_);
        residuals_.at(i).curWrap_.SetElementVector(&curLinState_);
        residuals_.at(i).noiWrap_.SetElementVector(&noise_);
        residuals_.at(i).res_->Eval(&residuals_.at(i).inn_,
                                    residuals_.at(i).preWrap_,
                                    residuals_.at(i).curWrap_,
                                    residuals_.at(i).noiWrap_);
        residuals_.at(i).innRef_.BoxMinus(residuals_.at(i).inn_,
                                          y.block(count, 0, residuals_.at(i).innDim_, 1));
        residuals_.at(i).res_->JacPre(residuals_.at(i).jacPre_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noiWrap_);
        residuals_.at(i).res_->JacCur(residuals_.at(i).jacCur_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noiWrap_);
        residuals_.at(i).res_->JacNoi(residuals_.at(i).jacNoi_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noiWrap_);
        residuals_.at(i).preWrap_.EmbedJacobian(JacPre,
                                                residuals_.at(i).jacPre_,
                                                count);
        residuals_.at(i).curWrap_.EmbedJacobian(JacCur,
                                                residuals_.at(i).jacCur_,
                                                count);
        residuals_.at(i).noiWrap_.EmbedJacobian(JacNoi,
                                                residuals_.at(i).jacNoi_,
                                                count);
        count += residuals_.at(i).innDim_;
      }
    }
    LOG(INFO) << "Innovation:\t" << y.transpose();

    // Compute weighting // TODO: make more efficient
    Winv = (JacNoi*R*JacNoi.transpose()).llt().solve(GIF::MatX::Identity(innDim,innDim));

    // Compute Kalman Update // TODO: make more efficient and numerically stable
    MatX D = cov_.inverse() + JacPre.transpose() * Winv * JacPre;
    MatX S = JacCur.transpose() * (Winv - Winv * JacPre * D.inverse() * JacPre.transpose() * Winv);
    cov_ = (S * JacCur).inverse();
    VecX dx = cov_ * S * y;

    // Apply Kalman Update
    curLinState_.BoxPlus(dx, &state_);
    time_ = t;
    LOG(INFO) << "state after Update:";
    LOG(INFO) << state_.Print();

    // Post Processing
    for (int i = 0; i < residuals_.size(); i++) {
      if (hasMeas.at(i)) {
        PostProcess(i);

        // Remove processed measurement
        if(residuals_.at(i).mt_.GetFirstTime() == t){
          residuals_.at(i).mt_.RemoveProcessedFirst();
        } else {
          LOG(WARNING) << "Bad timing";
        }
      }
    }
  }

  void TestJacs(const double delta, const double th, int i){
    residuals_.at(i).res_->TestJacs(delta, th);
  }

  void TestJacs(const double delta, const double th){
    for(int i=0;i<residuals_.size();i++){
      TestJacs(delta, th, i);
    }
  }

  virtual void PreProcess(int residual_id){};
  virtual void PostProcess(int residual_id){};

  ElementVector& GetState(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing state before initialization";
    return state_;
  }

  GIF::MatX& GetCovariance(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing cov before initialization";
    return cov_;
  }

  TimePoint& GetTime(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing time before initialization";
    return time_;
  }

  bool IsInitialized(){
    return is_initialized_;
  }

 protected:
  std::shared_ptr<ElementVectorDefinition> stateDefinition_; // Must come before state
  std::shared_ptr<ElementVectorDefinition> noiseDefinition_;
  std::vector<ResidualStruct> residuals_;
  TimePoint time_;
  TimePoint startTime_;
  ElementVector state_;
  ElementVector noise_;
  ElementVector curLinState_;
  MatX cov_;
  bool is_initialized_;

};

}

#endif /* GIF_FILTER_HPP_ */
