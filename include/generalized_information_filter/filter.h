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
  ResidualStruct(const BinaryResidualBase::Ptr& res,
                 const ElementVectorDefinition::Ptr& stateDefinition,
                 const ElementVectorDefinition::Ptr& noiseDefinition,
                 const Duration& maxWaitTime,
                 const Duration& minWaitTime):
                     mt_(!res->isUnary_, maxWaitTime, minWaitTime),
                     preWrap_(res->PreDefinition(), stateDefinition),
                     curWrap_(res->CurDefinition(), stateDefinition),
                     noiWrap_(res->NoiDefinition(), noiseDefinition),
                     inn_(res->InnDefinition()),
                     innRef_(res->InnDefinition()){
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
    isActive_ = false;
  }
  ~ResidualStruct() {
  }

  BinaryResidualBase::Ptr res_;
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
  bool isActive_;
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
    LOG(INFO) << "Initializing state at t = " << Print(t) << std::endl;
    startTime_ = t;
    time_ = t;
    if(initState != nullptr){
      state_ = *initState;
    } else {
      state_.SetIdentity();
    }
    inf_.setIdentity();
    is_initialized_ = true;
  }

  void Construct(){
    state_.Construct();
    state_.SetIdentity();
    curLinState_.Construct();
    noise_.Construct();
    noise_.SetIdentity();
    inf_.resize(stateDefinition_->GetDim(), stateDefinition_->GetDim());
  }

  int AddResidual(const BinaryResidualBase::Ptr& res,
                  const Duration& maxWaitTime,
                  const Duration& minWaitTime) {
    residuals_.emplace_back(res, stateDefinition_, noiseDefinition_, maxWaitTime, minWaitTime);
    Construct();
    return residuals_.size() - 1;
  }

  ElementVectorDefinition::Ptr StateDefinition() const {
    return stateDefinition_;
  }

  ElementVectorDefinition::Ptr NoiseDefinition() const {
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

  TimePoint GetMaxMinMeasTime() const {
    TimePoint maxMinMeasTime = TimePoint::min();
    bool check = true;
    for (int i = 0; i < residuals_.size(); i++) {
      if(!residuals_.at(i).res_->isUnary_){
        check = check & (residuals_.at(i).mt_.GetLastProcessedTime() != TimePoint::min());
        maxMinMeasTime = std::max(maxMinMeasTime, residuals_.at(i).mt_.GetLastProcessedTime());
      }
    }
    return check ? maxMinMeasTime : TimePoint::min();
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

  void AddMeasurement(const int i, const ElementVectorBase::CPtr& meas,
                      const TimePoint& t) {
    if(t <= time_){
      LOG(WARNING) << "Adding measurements before current time (will be discarded)" << std::endl;
      return;
    }
    LOG(INFO) << "Adding measurement with ID " << i << " at  t = " << Print(t) << std::endl;
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
    LOG(INFO) << out.str();
    for (int i = 0; i < residuals_.size(); i++) {
      LOG(INFO) << residuals_.at(i).mt_.Print(start, startOffset, resolution);
    }
  }

  void Update() {
    // Initialize if possible
    if(!is_initialized_ && GetMaxMinMeasTime() != TimePoint::min()){
      Init(GetMaxMinMeasTime());
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
      LOG(INFO) << "stateTime:\t" << Print(time_);
      TimePoint currentTime = GetCurrentTimeFromMeasurements();
      LOG(INFO) << "currentTime:\t" << Print(currentTime);
      TimePoint maxUpdateTime = GetMaxUpdateTime(currentTime);
      LOG(INFO) << "maxUpdateTime:\t" << Print(maxUpdateTime);
      std::set<TimePoint> times;
      GetMeasurementTimeList(times, maxUpdateTime, false);
      std::ostringstream out;
      out << "updateTimes:\t";
      for (const auto& t : times) {
        out << Print(t) << "\t";
      }
      LOG(INFO) << out.str();
      SplitAndMergeMeasurements(times);
      PrintMeasurementTimelines(time_, 20, 0.001);

      // Carry out updates
      for (const auto& t : times) {
        MakeUpdateStep(t);
      }
    }
  }

  void MakeUpdateStep(const TimePoint& t) {
    // Compute linearisation point // TODO: allow generic initialization
    curLinState_ = state_;

    // Check available measurements and prepare residuals
    int innDim = 0;
    ElementVectorBase::CPtr meas;
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).isActive_ = residuals_.at(i).mt_.GetMeasurement(t, meas);
      if (residuals_.at(i).isActive_) {
        residuals_.at(i).res_->SetDt(toSec(t-time_));
        residuals_.at(i).res_->SetMeas(meas);
        innDim += residuals_.at(i).res_->InnDefinition()->GetDim();
      }
    }
    PreProcess();

    // Temporaries
    VecX y(innDim);
    y.setZero();
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
      if (residuals_.at(i).isActive_) {
        ResidualStruct rs = residuals_.at(i);
        rs.preWrap_.SetElementVector(&state_);
        rs.curWrap_.SetElementVector(&curLinState_);
        rs.noiWrap_.SetElementVector(&noise_);
        rs.res_->Eval(&rs.inn_, rs.preWrap_,
                                rs.curWrap_,
                                rs.noiWrap_);
        rs.inn_.BoxMinus(rs.innRef_, y.block(count, 0, rs.innDim_, 1));
        LOG_IF(ERROR,y.hasNaN()) << "Residual " << rs.res_->name_
                                 << " contains NaN!\n" << y.transpose();
        rs.res_->JacPre(rs.jacPre_, rs.preWrap_,
                                    rs.curWrap_,
                                    rs.noiWrap_);
        LOG_IF(ERROR,rs.jacPre_.hasNaN()) << "jacPre " << rs.res_->name_ << " contains NaN!";
        rs.res_->JacCur(rs.jacCur_, rs.preWrap_,
                                    rs.curWrap_,
                                    rs.noiWrap_);
        LOG_IF(ERROR,rs.jacCur_.hasNaN()) << "jacCur " << rs.res_->name_ << " contains NaN!";
        rs.res_->JacNoi(rs.jacNoi_, rs.preWrap_,
                                    rs.curWrap_,
                                    rs.noiWrap_);
        LOG_IF(ERROR,rs.jacNoi_.hasNaN()) << "jacNoi " << rs.res_->name_ << " contains NaN!";
        for(int j=0;j<rs.res_->InnDefinition()->GetNumElements();j++){
          const ElementDescriptionBase::CPtr& description =
              rs.res_->InnDefinition()->GetElementDescription(j);
          if(!description->IsVectorSpace()){
            rs.jacPre_.block(rs.res_->InnDefinition()->GetStart(j),0,
                description->GetDim(),rs.preWrap_.GetDim()) =
                    rs.inn_.GetElement(j)->BoxminusJacInp(*rs.innRef_.GetElement(j)) *
                    rs.jacPre_.block(rs.res_->InnDefinition()->GetStart(j),0,
                                    description->GetDim(),rs.preWrap_.GetDim());
          }
        }

        rs.preWrap_.EmbedJacobian(JacPre, rs.jacPre_, count);
        rs.curWrap_.EmbedJacobian(JacCur, rs.jacCur_, count);
        rs.noiWrap_.EmbedJacobian(JacNoi, rs.jacNoi_, count);

        for(int j=0;j<rs.res_->NoiDefinition()->GetNumElements();j++){
          const ElementDescriptionBase::CPtr& description =
              rs.res_->NoiDefinition()->GetElementDescription(j);
          const int noiDim = description->GetDim();
          const int start1 = rs.res_->NoiDefinition()->GetStart(j);
          const int start2 = NoiseDefinition()->GetStart(NoiseDefinition()->FindName(rs.res_->NoiDefinition()->GetName(j)));
          R.block(start2,start2,noiDim,noiDim) = rs.res_->GetNoiseCovariance().block(start1,start1,noiDim,noiDim); // TODO: speed up and clean up
        }
        count += rs.innDim_;
      }
    }
    LOG(INFO) << "Innovation:\t" << y.transpose();

    // Compute weighting // TODO: make more efficient, exploit sparsity
    Eigen::LDLT<MatX> W_LDLT(JacNoi*R*JacNoi.transpose());
    LOG_IF(ERROR,W_LDLT.info() != Eigen::Success) << "Computation of Winv failed";
    Winv = W_LDLT.solve(GIF::MatX::Identity(innDim,innDim));

    // Compute Kalman Update
    MatX D = inf_ + JacPre.transpose() * Winv * JacPre;
    MatX S = JacCur.transpose() * (Winv - Winv * JacPre * D.inverse() * JacPre.transpose() * Winv);
    inf_ = S * JacCur;
    Eigen::LDLT<MatX> I_LDLT(inf_);
    LOG_IF(ERROR,I_LDLT.info() != Eigen::Success) << "Computation of Iinv failed";
    VecX dx = -I_LDLT.solve(S * y);

    // Apply Kalman Update
    curLinState_.BoxPlus(dx, &state_);
    time_ = t;
    LOG(INFO) << "state after Update:";
    LOG(INFO) << state_.Print();

    // Post Processing
    PostProcess();
    for (int i = 0; i < residuals_.size(); i++) {
      if (residuals_.at(i).isActive_) {
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

  virtual void PreProcess(){};
  virtual void PostProcess(){};

  ElementVector& GetState(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing state before initialization";
    return state_;
  }

  ElementVector& GetLinState(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing state before initialization";
    return curLinState_;
  }

  MatX GetCovariance(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing cov before initialization";
    const int stateDim = stateDefinition_->GetDim();
    return (inf_).ldlt().solve(GIF::MatX::Identity(stateDim,stateDim));
  }

  void SetCovariance(MatX cov){
    LOG_IF(ERROR,!is_initialized_) << "Accessing cov before initialization";
    const int stateDim = stateDefinition_->GetDim();
    inf_ = (cov).ldlt().solve(GIF::MatX::Identity(stateDim,stateDim));
  }

  inline MatRefX GetNoiseInfBlock(int i){
    const int dim = StateDefinition()->GetElementDescription(i)->GetDim();
    const int start = StateDefinition()->GetStart(i);
    return inf_.block(start,start,dim,dim);
  }

  inline MatRefX GetNoiseInfBlock(const std::string& str){
    const int outer_index = StateDefinition()->FindName(str);
    LOG_IF(FATAL, outer_index == -1) << "No element with name "
                                     << str << " for GetNoiseInfBlock!";
    return GetNoiseInfBlock(outer_index);
  }

  TimePoint& GetTime(){
    LOG_IF(ERROR,!is_initialized_) << "Accessing time before initialization";
    return time_;
  }

  bool IsInitialized(){
    return is_initialized_;
  }

  std::string PrintConnectivity(){
    std::ostringstream out;
    int resNamePadding = 10;

    std::vector<long> lengthState;
    for (int i = 0; i < state_.GetNumElement(); i++) {
      lengthState.push_back(state_.GetName(i).length());
      out << state_.GetName(i) << " ";
    }
    out << std::string(resNamePadding, '-') << " ";
    for (int i = 0; i < state_.GetNumElement(); i++) {
      out << state_.GetName(i) << " ";
    }
    out << std::endl;

    for (int i = 0; i < residuals_.size(); i++) {
      ResidualStruct rs = residuals_.at(i);
      for (int i = 0; i < state_.GetNumElement(); i++) {
        const char c = rs.preWrap_.FindName(state_.GetName(i)) == -1 ? ' ' : '+';
        out << std::string(lengthState.at(i), c);
        out << " ";
      }
      const std::string resName = rs.res_->name_.substr(0,resNamePadding);
      const int frontPadding = (resNamePadding-resName.length())/2;
      const int endPadding = resNamePadding-resName.length()-frontPadding+1;
      out << std::string(frontPadding, ' ') << resName << std::string(endPadding, ' ');
      for (int i = 0; i < state_.GetNumElement(); i++) {
        const char c = rs.curWrap_.FindName(state_.GetName(i)) == -1 ? ' ' : '+';
        out << std::string(lengthState.at(i), c);
        out << " ";
      }
      out << std::endl;
    }
    return out.str();
  }

 protected:
  ElementVectorDefinition::Ptr stateDefinition_; // Must come before state
  ElementVectorDefinition::Ptr noiseDefinition_;
  std::vector<ResidualStruct> residuals_;
  TimePoint time_;
  TimePoint startTime_;
  ElementVector state_;
  ElementVector noise_;
  ElementVector curLinState_;
  MatX inf_;
  bool is_initialized_;

};

}

#endif /* GIF_FILTER_HPP_ */
