#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/measurement.h"

namespace GIF {

class ResidualStruct {
 public:
  ResidualStruct(const std::shared_ptr<BinaryResidualBase>& res,
                 const std::shared_ptr<ElementVectorDefinition>& stateDefinition) {
    stateDefinition->Extend(res->preDefinition());
    stateDefinition->Extend(res->curDefinition());
    res_ = res;
    mt_.reset(new MeasurementTimeline(!res->isUnary_));
    preWrap_.reset(new ElementVectorWrapper(res->preDefinition(), stateDefinition));
    curWrap_.reset(new ElementVectorWrapper(res->curDefinition(), stateDefinition));
    inn_.reset(new ElementVector(res->innDefinition()));
    innRef_.reset(new ElementVector(res->innDefinition()));
    innRef_->SetIdentity();
    noi_.reset(new ElementVector(res->noiDefinition()));
    noi_->SetIdentity();
    innDim_ = res->innDefinition()->GetStateDimension();
    jacPre_.resize(innDim_, res->preDefinition()->GetStateDimension());
    jacCur_.resize(innDim_, res->curDefinition()->GetStateDimension());
    jacNoi_.resize(innDim_, res->noiDefinition()->GetStateDimension());
  }
  ;
  ~ResidualStruct() {
  }
  ;

  std::shared_ptr<BinaryResidualBase> res_;
  std::shared_ptr<MeasurementTimeline> mt_;
  std::shared_ptr<ElementVectorWrapper> preWrap_;
  std::shared_ptr<ElementVectorWrapper> curWrap_;
  std::shared_ptr<ElementVector> inn_;
  std::shared_ptr<ElementVector> innRef_;
  std::shared_ptr<ElementVector> noi_;
  MatX jacPre_;
  MatX jacCur_;
  MatX jacNoi_;
  int innDim_;
};

class Filter {
 public:
  Filter()
      : stateDefinition_(new ElementVectorDefinition()) {
    init();
  }
  ;
  virtual ~Filter() {
  }
  ;

  void init(const TimePoint& t = TimePoint::min()) {  // TODO: pass optional State
    startTime_ = t;
    time_ = t;
    state_.reset(new ElementVector(stateDefinition_));
    state_->SetIdentity();
    curLinState_.reset(new ElementVector(stateDefinition_));
    cov_.resize(stateDefinition_->GetStateDimension(), stateDefinition_->GetStateDimension());
    cov_.setIdentity();
  }

  int addRes(const std::shared_ptr<BinaryResidualBase>& res,
             const std::string& name = "") {
    residuals_.emplace_back(res, stateDefinition_);
    return residuals_.size() - 1;
  }
  void evalRes(const std::shared_ptr<const ElementVectorBase>& pre,
               const std::shared_ptr<const ElementVectorBase>& cur) {
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).preWrap_->SetVectorElement(pre);
      residuals_.at(i).curWrap_->SetVectorElement(cur);
      residuals_.at(i).res_->eval(residuals_.at(i).inn_,
                                          residuals_.at(i).preWrap_,
                                          residuals_.at(i).curWrap_,
                                          residuals_.at(i).noi_);
      residuals_.at(i).inn_->Print();
    }
  }
  std::shared_ptr<ElementVectorDefinition> stateDefinition() const {
    return stateDefinition_;
  }
  TimePoint getCurrentTimeFromMeasurements() const {
    TimePoint currentTime = TimePoint::min();
    for (int i = 0; i < residuals_.size(); i++) {
      currentTime = std::max(currentTime, residuals_.at(i).mt_->getLastTime());
    }
    return currentTime;
  }
  TimePoint getMaxUpdateTime(const TimePoint& currentTime) const {
    TimePoint maxUpdateTime = currentTime;
    for (int i = 0; i < residuals_.size(); i++) {
      maxUpdateTime = std::min(
          maxUpdateTime,
          residuals_.at(i).mt_->getMaximalUpdateTime(currentTime));
    }
    return maxUpdateTime;
  }
  void getMeasurementTimeList(std::set<TimePoint>& times,
                              const TimePoint& maxUpdateTime,
                              const bool includeMax = false) const {
    for (int i = 0; i < residuals_.size(); i++) {
      if (!residuals_.at(i).res_->isMergeable_) {
        // Add all non-mergeable measurement times
        residuals_.at(i).mt_->getAllInRange(times, time_, maxUpdateTime);
      } else if (!residuals_.at(i).res_->isSplitable_
          && residuals_.at(i).res_->isUnary_) {
        // For the special case of unary and mergeable residuals add the last measurement time
        residuals_.at(i).mt_->getLastInRange(times, time_, maxUpdateTime);
      }
    }
    if (includeMax) {
      times.insert(maxUpdateTime);
    }
  }
  void splitAndMergeMeasurements(const std::set<TimePoint>& times) {
    for (int i = 0; i < residuals_.size(); i++) {
      if (residuals_.at(i).res_->isSplitable_
          && !residuals_.at(i).res_->isUnary_) {
        // First insert all (splitable + !unary)
        residuals_.at(i).mt_->split(times, residuals_.at(i).res_);
      }
      if (residuals_.at(i).res_->isMergeable_) {
        // Then remove all undesired (mergeable)
        residuals_.at(i).mt_->mergeUndesired(times, residuals_.at(i).res_);
      }
    }
  }
  void addMeas(const int i, const std::shared_ptr<const ElementVectorBase>& meas,
               const TimePoint& t) {
    residuals_.at(i).mt_->addMeas(meas, t);
  }
  void printMeasurementTimelines(const TimePoint& start = TimePoint::min(),
                                 int startOffset = 0,
                                 double resolution = 0.01) {
    for (int i = 0; i < startOffset; i++) {
      std::cout << " ";
    }
    std::cout << "|" << std::endl;
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).mt_->print(start, startOffset, resolution);
    }
  }
  void update() {
    // Remove outdated
    for (int i = 0; i < residuals_.size(); i++) {
      residuals_.at(i).mt_->removeOutdated(time_);
    }
    printMeasurementTimelines(time_, 20);
    std::cout << "stateTime:\t" << toSec(time_ - startTime_) << std::endl;
    TimePoint currentTime = getCurrentTimeFromMeasurements();
    std::cout << "currentTime:\t" << toSec(currentTime - startTime_)
        << std::endl;
    TimePoint maxUpdateTime = getMaxUpdateTime(currentTime);
    std::cout << "maxUpdateTime:\t" << toSec(maxUpdateTime - startTime_)
        << std::endl;
    std::set<TimePoint> times;
    getMeasurementTimeList(times, maxUpdateTime, false);
    std::cout << "updateTimes:\t";
    for (auto t : times) {
      std::cout << toSec(t - startTime_) << "\t";
    }
    std::cout << std::endl;
    splitAndMergeMeasurements(times);
    printMeasurementTimelines(time_, 20);

    for (auto t : times) {
      makeUpdateStep(t);
    }
  }

  void makeUpdateStep(const TimePoint& t) {
    // Compute linearisation point
    curLinState_ = state_;

    // Eval residual and Jacobians
    int innDim = 0;
    std::vector<bool> hasMeas(residuals_.size(), false);
    std::shared_ptr<const ElementVectorBase> meas;
    for (int i = 0; i < residuals_.size(); i++) {
      hasMeas.at(i) = residuals_.at(i).mt_->getMeas(t, meas);
      if (hasMeas.at(i)) {
        innDim += residuals_.at(i).res_->innDefinition()->GetStateDimension();
      }
    }
    VecX y(innDim);
    MatX jacPre(innDim, stateDefinition_->GetStateDimension());
    MatX jacCur(innDim, stateDefinition_->GetStateDimension());
    jacPre.setZero();
    jacCur.setZero();
    MatX Winv(innDim, innDim);
    Winv.setZero();
    int count = 0;
    for (int i = 0; i < residuals_.size(); i++) {
      if (hasMeas.at(i)) {
        residuals_.at(i).mt_->getMeas(t, meas);
        residuals_.at(i).res_->setMeas(meas);
        residuals_.at(i).preWrap_->SetVectorElement(state_);
        residuals_.at(i).curWrap_->SetVectorElement(curLinState_);
//        residuals_.at(i).res_->testJacs(residuals_.at(i).preWrap_,residuals_.at(i).curWrap_,residuals_.at(i).noi_);
        residuals_.at(i).res_->eval(residuals_.at(i).inn_,
                                            residuals_.at(i).preWrap_,
                                            residuals_.at(i).curWrap_,
                                            residuals_.at(i).noi_);
        residuals_.at(i).innRef_->BoxMinus(
            residuals_.at(i).inn_,
            y.block(count, 0, residuals_.at(i).innDim_, 1));

        // Compute Jacobians
        residuals_.at(i).res_->jacPre(residuals_.at(i).jacPre_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).res_->jacCur(residuals_.at(i).jacCur_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).res_->jacNoi(residuals_.at(i).jacNoi_,
                                      residuals_.at(i).preWrap_,
                                      residuals_.at(i).curWrap_,
                                      residuals_.at(i).noi_);
        residuals_.at(i).preWrap_->EmbedJacobian(jacPre,
                                                residuals_.at(i).jacPre_,
                                                count);
        residuals_.at(i).curWrap_->EmbedJacobian(jacCur,
                                                residuals_.at(i).jacCur_,
                                                count);
        Winv.block(count, count, residuals_.at(i).innDim_,
                   residuals_.at(i).innDim_) = (residuals_.at(i).jacNoi_
            * residuals_.at(i).res_->getR()
            * residuals_.at(i).jacNoi_.transpose()).inverse();

        // Increment counter
        count += residuals_.at(i).innDim_;
      }
    }
    std::cout << "Innovation:\t" << y.transpose() << std::endl;
//    std::cout << jacPre << std::endl;
//    std::cout << jacCur << std::endl;
//    std::cout << Winv << std::endl;

    // Compute Kalman update // TODO: make more efficient
    MatX D = cov_.inverse() + jacPre.transpose() * Winv * jacPre;
    MatX S = jacCur.transpose()
        * (Winv - Winv * jacPre * D.inverse() * jacPre.transpose() * Winv);
    cov_ = (S * jacCur).inverse();
    VecX dx = cov_ * S * y;

    // Apply Kalman update
    curLinState_->BoxPlus(dx, state_);
    std::cout << "state after update:" << std::endl;
    state_->Print();

    // Remove measurements and update timepoint
    time_ = t;
  }

 protected:
  TimePoint time_;
  TimePoint startTime_;
  std::shared_ptr<ElementVector> state_;
  std::shared_ptr<ElementVector> curLinState_;
  MatX cov_;

  std::shared_ptr<ElementVectorDefinition> stateDefinition_;
  std::vector<ResidualStruct> residuals_;
};

}

#endif /* GIF_FILTER_HPP_ */
