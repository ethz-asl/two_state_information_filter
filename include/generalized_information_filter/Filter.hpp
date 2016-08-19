/*
 * Filter.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include "generalized_information_filter/BinaryResidual.hpp"
#include "generalized_information_filter/common.hpp"

namespace GIF{

class Filter{
 public:
  Filter(): stateDefinition_(new StateDefinition()){
    init();
  };
  virtual ~Filter(){};

  void init(const TimePoint& t = TimePoint::min()){ // TODO: pass optional State
    startTime_ = t;
    time_ = t;
    state_.reset(new State(stateDefinition_));
    posLinState_.reset(new State(stateDefinition_));
    cov_.setIdentity();
  }

  void addRes(const std::shared_ptr<BinaryResidualBase>& r, const std::string& name = ""){
    binaryResiduals_.push_back(r);
    stateDefinition_->extend(r->preDefinition());
    stateDefinition_->extend(r->posDefinition());
    binaryMeasurementTimelines_.emplace_back(new MeasurementTimeline(!r->isUnary_));
    binaryWrappersPre_.emplace_back(new StateWrapper(r->preDefinition(),stateDefinition_));
    binaryWrappersPos_.emplace_back(new StateWrapper(r->posDefinition(),stateDefinition_));
  }
  void evalRes(const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos){
    for(int i=0;i<binaryResiduals_.size();i++){
      std::shared_ptr<State> inn(new State(binaryResiduals_.at(i)->resDefinition()));
      std::shared_ptr<State> noi(new State(binaryResiduals_.at(i)->noiDefinition()));
      noi->setIdentity();
      binaryWrappersPre_.at(i)->setState(pre);
      binaryWrappersPos_.at(i)->setState(pos);
      binaryResiduals_.at(i)->evalResidual(inn,binaryWrappersPre_.at(i),binaryWrappersPos_.at(i),noi);
      inn->print();
    }
  }
  std::shared_ptr<StateDefinition> stateDefinition() const{
    return stateDefinition_;
  }
  TimePoint getCurrentTimeFromMeasurements() const{
    TimePoint currentTime = TimePoint::min();
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      currentTime = std::max(currentTime,binaryMeasurementTimelines_.at(i)->getLastTime());
    }
    return currentTime;
  }
  TimePoint getMaxUpdateTime(const TimePoint& currentTime) const{
    TimePoint maxUpdateTime = currentTime;
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      maxUpdateTime = std::min(maxUpdateTime,binaryMeasurementTimelines_.at(i)->getMaximalUpdateTime(currentTime));
    }
    return maxUpdateTime;
  }
  void getMeasurementTimeList(std::set<TimePoint>& times, const TimePoint& maxUpdateTime, const bool includeMax = false) const{
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      if(!binaryResiduals_.at(i)->isMergeable_){
        // Add all non-mergeable measurement times
        binaryMeasurementTimelines_.at(i)->addAllInRange(times,time_,maxUpdateTime);
      } else if(!binaryResiduals_.at(i)->isSplitable_ && binaryResiduals_.at(i)->isUnary_){
        // For the special case of unary and mergeable residuals add the last measurement time
        binaryMeasurementTimelines_.at(i)->addLastInRange(times,time_,maxUpdateTime);
      }
    }
    if(includeMax){
      times.insert(maxUpdateTime);
    }
  }
  void splitAndMergeMeasurements(const std::set<TimePoint>& times){
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      if(binaryResiduals_.at(i)->isSplitable_ && !binaryResiduals_.at(i)->isUnary_){
        // First insert all (splitable + !unary)
        binaryMeasurementTimelines_.at(i)->split(times,binaryResiduals_.at(i));
      }
      if(binaryResiduals_.at(i)->isMergeable_){
        // Then remove all undesired (mergeable)
        binaryMeasurementTimelines_.at(i)->mergeUndesired(times,binaryResiduals_.at(i));
      }
    }
  }
  void addMeas(const int i, const std::shared_ptr<const MeasurementBase>& meas, const TimePoint& t){
    binaryMeasurementTimelines_.at(i)->addMeas(meas,t);
  }
  void printMeasurementTimelines(const TimePoint& start = TimePoint::min()){
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      std::cout << "Timeline " << i << std::endl;
      binaryMeasurementTimelines_.at(i)->print(start);
    }
  }
  void update(){
    printMeasurementTimelines(startTime_);
    std::cout << "stateTime " << toSec(time_-startTime_) << std::endl;
    TimePoint currentTime = getCurrentTimeFromMeasurements();
    std::cout << "currentTime " << toSec(currentTime-startTime_) << std::endl;
    TimePoint maxUpdateTime = getMaxUpdateTime(currentTime);
    std::cout << "maxUpdateTime " << toSec(maxUpdateTime-startTime_) << std::endl;
    std::set<TimePoint> times;
    getMeasurementTimeList(times,maxUpdateTime,false);
    std::cout << "updateTimes ";
    for(auto t : times){
      std::cout << toSec(t-startTime_) << "\t";
    }
    std::cout << std::endl;
    splitAndMergeMeasurements(times);
    printMeasurementTimelines(startTime_);

    for(auto t : times){
      makeUpdateStep(t);
    }
  }

  void makeUpdateStep(const TimePoint& t){
    // Compute linearisation point
    posLinState_ = state_;

    // Eval residual and Jacobians
    std::vector<bool> hasMeas(binaryResiduals_.size(),false);
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      std::shared_ptr<const MeasurementBase> meas;
      hasMeas.at(i) = binaryMeasurementTimelines_.at(i)->getMeas(t,meas);
      if(hasMeas.at(i)){
        binaryResiduals_.at(i)->setMeas(meas);
        std::shared_ptr<State> inn(new State(binaryResiduals_.at(i)->resDefinition()));
        std::shared_ptr<State> noi(new State(binaryResiduals_.at(i)->noiDefinition()));
        noi->setIdentity();
        binaryWrappersPre_.at(i)->setState(state_);
        binaryWrappersPos_.at(i)->setState(posLinState_);
        binaryResiduals_.at(i)->evalResidual(inn,binaryWrappersPre_.at(i),binaryWrappersPos_.at(i),noi);
        inn->print();
      }
    }


    // Compute Kalman update

    // Apply Kalman update

    // Remove measurements and update timepoint
    time_ = t;
  }

 protected:
  TimePoint time_;
  TimePoint startTime_;
  std::shared_ptr<State> state_;
  std::shared_ptr<State> posLinState_;
  MXD cov_;

  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
  std::vector<std::shared_ptr<MeasurementTimeline>> binaryMeasurementTimelines_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPre_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPos_;
};

}

#endif /* GIF_FILTER_HPP_ */
