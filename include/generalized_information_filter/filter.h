#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/measurement.h"

namespace GIF {

/*! \brief Residual Struct
 *         Contains various object and temporaries associated with a specific residual.
 */
class ResidualStruct {
 public:
  ResidualStruct(const std::shared_ptr<BinaryResidualBase>& res,
                 const std::shared_ptr<ElementVectorDefinition>& stateDefinition,
                 const Duration& maxWaitTime,
                 const Duration& minWaitTime):
                     mt_(!res->isUnary_, maxWaitTime, minWaitTime),
                     preWrap_(res->PreDefinition(), stateDefinition),
                     curWrap_(res->CurDefinition(), stateDefinition),
                     inn_(res->InnDefinition()),
                     innRef_(res->InnDefinition()),
                     noi_(res->NoiDefinition()){
    stateDefinition->ExtendWithOtherElementVectorDefinition(*res->PreDefinition());
    stateDefinition->ExtendWithOtherElementVectorDefinition(*res->CurDefinition());
    res_ = res;
    preWrap_.ComputeMap();
    curWrap_.ComputeMap();
    innRef_.SetIdentity();
    noi_.SetIdentity();
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
  ElementVector inn_;
  ElementVector innRef_;
  ElementVector noi_;
  MatX jacPre_;
  MatX jacCur_;
  MatX jacNoi_;
  int innDim_;
};

/*! \brief Filter
 *         Handles residuals. Builds state defintion. Contains measurements timelines. Implements
 *         timing logic.
 */
class Filter {
 public:
  Filter(): stateDefinition_(new ElementVectorDefinition()),
            state_(stateDefinition_), curLinState_(stateDefinition_) {
    Init();
    max_wait_time_default_ = fromSec(0.1);
    min_wait_time_default_ = fromSec(0.0);
  }

  virtual ~Filter() {
  }

  void Init(const TimePoint& t = TimePoint::min()) {  // TODO: pass optional State
    startTime_ = t;
    time_ = t;
    state_.Construct();
    state_.SetIdentity();
    curLinState_.Construct();
    cov_.resize(stateDefinition_->GetDim(), stateDefinition_->GetDim());
    cov_.setIdentity();
  }

  int AddResidual(const std::shared_ptr<BinaryResidualBase>& res, const std::string& name = "") {
    residuals_.emplace_back(res, stateDefinition_,max_wait_time_default_,min_wait_time_default_);
    return residuals_.size() - 1;
  }

  void EvalResidual(const ElementVectorBase* pre, const ElementVectorBase* cur) {
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).preWrap_.SetElementVector(pre);
      residuals_.at(i).curWrap_.SetElementVector(cur);
      residuals_.at(i).res_->Eval(&residuals_.at(i).inn_,
                                  residuals_.at(i).preWrap_,
                                  residuals_.at(i).curWrap_,
                                  residuals_.at(i).noi_);
      residuals_.at(i).inn_.Print();
    }
  }

  std::shared_ptr<ElementVectorDefinition> StateDefinition() const {
    return stateDefinition_;
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

  void GetMeasurementTimeList(std::set<TimePoint>& times,
                              const TimePoint& maxUpdateTime,
                              const bool includeMax = false) const {
    for (int i = 0; i < residuals_.size(); i++) {
      if (!residuals_.at(i).res_->isMergeable_) {
        // Add all non-mergeable measurement times
        residuals_.at(i).mt_.GetAllInRange(times, time_, maxUpdateTime);
      } else if (!residuals_.at(i).res_->isSplitable_ && residuals_.at(i).res_->isUnary_) {
        // For the special case of unary and mergeable residuals add the last measurement time
        residuals_.at(i).mt_.GetLastInRange(times, time_, maxUpdateTime);
      }
    }
    if (includeMax) {
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
    residuals_.at(i).mt_.AddMeasurement(meas, t);
  }

  void PrintMeasurementTimelines(const TimePoint& start, int startOffset, double resolution) {
    for (int i = 0; i < startOffset; i++) {
      std::cout << " ";
    }
    std::cout << "|" << std::endl;
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).mt_.Print(start, startOffset, resolution);
    }
  }

  void Update() {
    // Remove outdated
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).mt_.RemoveOutdated(time_);
    }
    PrintMeasurementTimelines(time_, 20, 0.01);
    std::cout << "stateTime:\t" << toSec(time_ - startTime_) << std::endl;
    TimePoint currentTime = GetCurrentTimeFromMeasurements();
    std::cout << "currentTime:\t" << toSec(currentTime - startTime_) << std::endl;
    TimePoint maxUpdateTime = GetMaxUpdateTime(currentTime);
    std::cout << "maxUpdateTime:\t" << toSec(maxUpdateTime - startTime_) << std::endl;
    std::set<TimePoint> times;
    GetMeasurementTimeList(times, maxUpdateTime, false);
    std::cout << "updateTimes:\t";
    for (const auto& t : times) {
      std::cout << toSec(t - startTime_) << "\t";
    }
    std::cout << std::endl;
    SplitAndMergeMeasurements(times);
    PrintMeasurementTimelines(time_, 20, 0.01);

    for (const auto& t : times) {
      MakeUpdateStep(t);
    }
  }

  void MakeUpdateStep(const TimePoint& t) {
    // Compute linearisation point
    curLinState_ = state_;

    // Eval residual and Jacobians
    int innDim = 0;
    std::vector<bool> hasMeas(residuals_.size(), false);
    std::shared_ptr<const ElementVectorBase> meas;
    for (int i = 0; i < residuals_.size(); i++) {
      hasMeas.at(i) = residuals_.at(i).mt_.GetMeasurement(t, meas);
      if (hasMeas.at(i)) {
        innDim += residuals_.at(i).res_->InnDefinition()->GetDim();
      }
    }
    VecX y(innDim);
    MatX JacPre(innDim, stateDefinition_->GetDim());
    MatX JacCur(innDim, stateDefinition_->GetDim());
    JacPre.setZero();
    JacCur.setZero();
    MatX Winv(innDim, innDim);
    Winv.setZero();
    int count = 0;
    for (int i = 0; i < residuals_.size(); i++) {
      if (hasMeas.at(i)) {
        residuals_.at(i).mt_.GetMeasurement(t, meas);
        residuals_.at(i).res_->SetMeas(meas);
        residuals_.at(i).preWrap_.SetElementVector(&state_);
        residuals_.at(i).curWrap_.SetElementVector(&curLinState_);
        residuals_.at(i).res_->Eval(&residuals_.at(i).inn_,
                                    residuals_.at(i).preWrap_,
                                    residuals_.at(i).curWrap_,
                                    residuals_.at(i).noi_);
        residuals_.at(i).innRef_.BoxMinus(residuals_.at(i).inn_,
                                          y.block(count, 0, residuals_.at(i).innDim_, 1));

        // Compute Jacobians
        residuals_.at(i).res_->JacPre(residuals_.at(i).jacPre_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).res_->JacCur(residuals_.at(i).jacCur_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).res_->JacNoi(residuals_.at(i).jacNoi_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).preWrap_.EmbedJacobian(JacPre,
                                                residuals_.at(i).jacPre_,
                                                count);
        residuals_.at(i).curWrap_.EmbedJacobian(JacCur,
                                                residuals_.at(i).jacCur_,
                                                count);
        Winv.block(count, count, residuals_.at(i).innDim_,
                   residuals_.at(i).innDim_) = (residuals_.at(i).jacNoi_
                                               * residuals_.at(i).res_->GetNoiseCovariance()
                                               * residuals_.at(i).jacNoi_.transpose()).inverse();

        // Increment counter
        count += residuals_.at(i).innDim_;
      }
    }
    std::cout << "Innovation:\t" << y.transpose() << std::endl;

    // Compute Kalman Update // TODO: make more efficient and numerically stable
    MatX D = cov_.inverse() + JacPre.transpose() * Winv * JacPre;
    MatX S = JacCur.transpose() * (Winv - Winv * JacPre * D.inverse() * JacPre.transpose() * Winv);
    cov_ = (S * JacCur).inverse();
    VecX dx = cov_ * S * y;

    // Apply Kalman Update
    curLinState_.BoxPlus(dx, &state_);
    std::cout << "state after Update:" << std::endl;
    state_.Print();

    // Remove measurements and Update timepoint
    time_ = t;
  }

 protected:
  std::shared_ptr<ElementVectorDefinition> stateDefinition_; // Must come before state
  std::vector<ResidualStruct> residuals_;
  TimePoint time_;
  TimePoint startTime_;
  ElementVector state_;
  ElementVector curLinState_;
  MatX cov_;
  Duration max_wait_time_default_;
  Duration min_wait_time_default_;

};

}

#endif /* GIF_FILTER_HPP_ */
