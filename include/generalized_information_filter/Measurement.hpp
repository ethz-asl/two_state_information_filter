/*
 * Measurement.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_

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
};


//template<typename Meas>
//class MeasurementTimeline{
// public:
//  typedef Meas mtMeas;
//  std::map<double,mtMeas> measMap_;
//  typename std::map<double,mtMeas>::iterator itMeas_;
//  double maxWaitTime_;
//  double minWaitTime_;
//  MeasurementTimeline(){
//    maxWaitTime_ = 0.1;
//    minWaitTime_ = 0.0;
//  };
//  virtual ~MeasurementTimeline(){};
//  void addMeas(const mtMeas& meas, const double& t){
//    measMap_[t] = meas;
//  }
//  void clear()
//  {
//    measMap_.clear();
//  }
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
