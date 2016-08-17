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
  Filter(): stateDefinition_(new StateDefinition()){};
  virtual ~Filter(){};

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
  bool getCurrentTimeFromMeasurements(double& currentTime) const{
    bool foundMeasurement;
    double lastMeasurementTime;
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
  double getMaxUpdateTime(const double& currentTime) const{
    double maxUpdateTime = currentTime;
    for(int i=0;i<binaryMeasurementTimelines_.size();i++){
      maxUpdateTime = std::min(maxUpdateTime,binaryMeasurementTimelines_[i]->getMaximalUpdateTime(currentTime));
    }
    return maxUpdateTime;
  }
  void getMeasurementTimeList(std::list<double>& times, const double& maxUpdateTime) const{
    // Compose list of times which need to be processed
    // First register all time between current state time and maxUpdateTime
    // Second check which times can be eliminated
  }
  void splitAndMergeMeasurements(const std::list<double>& times){
    // Change the measurement timelines by applying the split and merge functionality of the residuals
    // UnaryDiscrete: split only second, merge not possible
    // Differential: split leads to twice the same, merge as weighted mean
  }

 protected:
  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
  std::vector<std::shared_ptr<MeasurementTimeline>> binaryMeasurementTimelines_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPre_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPos_;
};

}

#endif /* GIF_FILTER_HPP_ */
