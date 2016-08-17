/*
 * Measurement.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_
#include <map>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

class MeasurementBase: public State{
 public:
  MeasurementBase(const std::shared_ptr<const StateDefinition>& def): State(def){
  }
  virtual ~MeasurementBase(){};
};

class MeasurementTimeline{
 public:
  MeasurementTimeline(const double& maxWaitTime = 0.1, const double& minWaitTime = 0.0){
    maxWaitTime_ = maxWaitTime;
    minWaitTime_ = minWaitTime;
    lastProcessedTime_ = 0.0;
    hasProcessedTime_ = false;
  };
  virtual ~MeasurementTimeline(){};
  void addMeas(const std::shared_ptr<const MeasurementBase>& meas, const double& t){
    if(hasProcessedTime_ && t<lastProcessedTime_){
      std::cout << "Error: adding measurements before last processed time (will be discarded)" << std::endl;
    } else {
      measMap_[t] = meas;
    }
  }
  void removeProcessedFirst(){
    assert(measMap_.size() > 0);
    lastProcessedTime_ = measMap_.begin()->first;
    measMap_.erase(measMap_.begin());
    hasProcessedTime_ = true;
  }
  void removeProcessedMeas(const double& t){
    assert(measMap_.count(t) > 0);
    measMap_.erase(t);
    lastProcessedTime_ = t;
    hasProcessedTime_ = true;
  }
  void clear(){
    measMap_.clear();
    hasProcessedTime_ = false;
  }
  bool getLastTime(double& lastTime) const{
    if(!measMap_.empty()){
      lastTime = measMap_.rbegin()->first;
      return true;
    } else if(hasProcessedTime_){
      lastTime = lastProcessedTime_;
      return true;
    } else {
      return false;
    }
  }
  double getMaximalUpdateTime(const double& currentTime) const{
    double maximalUpdateTime = currentTime-maxWaitTime_;
    if(!measMap_.empty()){
      maximalUpdateTime = std::max(maximalUpdateTime,measMap_.rbegin()->first+minWaitTime_);
    } else if(hasProcessedTime_){
      maximalUpdateTime = std::max(maximalUpdateTime,lastProcessedTime_+minWaitTime_);
    }
    return maximalUpdateTime;
  }
 protected:
  std::map<double,std::shared_ptr<const MeasurementBase>> measMap_;
  double maxWaitTime_;
  double minWaitTime_;
  double lastProcessedTime_;
  bool hasProcessedTime_;
};

//class MeasurementTimeline{
//  void clean(double t){
//    while(measMap_.size() > 1 && measMap_.begin()->first<=t){
//      measMap_.erase(measMap_.begin());
//    }
//  }
//  bool getNextTime(double actualTime, double& nextTime){
//    itMeas_ = measMap_.upper_bound(actualTime);
//    if(itMeas_!=measMap_.end()){
//      nextTime = itMeas_->first;
//      return true;
//    } else {
//      return false;
//    }
//  }
//  void waitTime(double actualTime, double& time){
//    double measurementTime = actualTime-maxWaitTime_;
//    if(!measMap_.empty() && measMap_.rbegin()->first + minWaitTime_ > measurementTime){
//      measurementTime = measMap_.rbegin()->first + minWaitTime_;
//    }
//    if(time > measurementTime){
//      time = measurementTime;
//    }
//  }
//  bool getLastTime(double& lastTime){
//    if(!measMap_.empty()){
//      lastTime = measMap_.rbegin()->first;
//      return true;
//    } else {
//      return false;
//    }
//  }
//  bool hasMeasurementAt(double t){
//    return measMap_.count(t)>0;
//  }
//};

}

#endif /* GIF_MEASUREMENT_HPP_ */
