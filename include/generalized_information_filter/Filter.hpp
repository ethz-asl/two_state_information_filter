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

  void init(const TimePoint& t = TimePoint::min()){
    time_ = t;
    state_.reset(new State(stateDefinition_));
    cov_.setIdentity();
  }

  void addRes(const std::shared_ptr<BinaryResidualBase>& r, const std::string& name = ""){
    binaryResiduals_.push_back(r);
    stateDefinition_->extend(r->preDefinition());
    stateDefinition_->extend(r->posDefinition());
    binaryMeasurementTimelines_.emplace_back(new MeasurementTimeline());
    binaryWrappersPre_.emplace_back(new StateWrapper(r->preDefinition(),stateDefinition_));
    binaryWrappersPos_.emplace_back(new StateWrapper(r->posDefinition(),stateDefinition_));
  }
  void evalRes(const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos){
    for(int i=0;i<binaryResiduals_.size();i++){
      std::shared_ptr<State> inn(new State(binaryResiduals_[i]->resDefinition()));
      std::shared_ptr<State> noi(new State(binaryResiduals_[i]->noiDefinition()));
      noi->setIdentity();
      binaryWrappersPre_[i]->setState(pre);
      binaryWrappersPos_[i]->setState(pos);
      binaryResiduals_[i]->evalResidual(inn,binaryWrappersPre_[i],binaryWrappersPos_[i],noi);
      inn->print();
    }
  }
  std::shared_ptr<StateDefinition> stateDefinition() const{
    return stateDefinition_;
  }
  bool getCurrentTimeFromMeasurements(TimePoint& currentTime) const{
    bool foundMeasurement;
    TimePoint lastMeasurementTime;
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      if(!foundMeasurement){
        foundMeasurement = binaryMeasurementTimelines_[i]->getLastTime(currentTime);
      } else {
        if(binaryMeasurementTimelines_[i]->getLastTime(lastMeasurementTime) && lastMeasurementTime > currentTime){
          currentTime = lastMeasurementTime;
        }
      }
    }
    return foundMeasurement;
  }
  TimePoint getMaxUpdateTime(const TimePoint& currentTime) const{
    TimePoint maxUpdateTime = currentTime;
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      maxUpdateTime = std::min(maxUpdateTime,binaryMeasurementTimelines_[i]->getMaximalUpdateTime(currentTime));
    }
    return maxUpdateTime;
  }
  void getMeasurementTimeList(std::set<TimePoint>& times, const TimePoint& maxUpdateTime, const bool includeMax = false) const{
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      if(!binaryResiduals_[i]->isMergeable_){
        binaryMeasurementTimelines_[i]->addAllInRange(times,time_,maxUpdateTime);
      } else if(!binaryResiduals_[i]->isSplitable_ && binaryResiduals_[i]->isUnary_){
        binaryMeasurementTimelines_[i]->addLastInRange(times,time_,maxUpdateTime);
      }
    }
    if(includeMax){
      times.insert(maxUpdateTime);
    }
  }
  void splitAndMergeMeasurements(const std::set<TimePoint>& times){
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      // First insert all (splitable + !unary)
      // Then remove all undisered (mergeable)
    }
  }

 protected:
  TimePoint time_;
  std::shared_ptr<State> state_;
  MXD cov_;

  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
  std::vector<std::shared_ptr<MeasurementTimeline>> binaryMeasurementTimelines_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPre_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPos_;
};

}

#endif /* GIF_FILTER_HPP_ */
