#include "generalized_information_filter/Measurement.hpp"
#include "generalized_information_filter/BinaryResidual.hpp"

namespace GIF{

MeasurementTimeline::MeasurementTimeline(const Duration& maxWaitTime, const Duration& minWaitTime){
  maxWaitTime_ = maxWaitTime;
  minWaitTime_ = minWaitTime;
  lastProcessedTime_ = TimePoint::min();
  hasProcessedTime_ = false;
};

MeasurementTimeline::~MeasurementTimeline(){};

void MeasurementTimeline::addMeas(const std::shared_ptr<const MeasurementBase>& meas, const TimePoint& t){
  if(hasProcessedTime_ && t<=lastProcessedTime_){
    std::cout << "Error: adding measurements before last processed time (will be discarded)" << std::endl;
  } else {
    std::pair<std::map<TimePoint,std::shared_ptr<const MeasurementBase>>::iterator,bool> ret;
    ret = measMap_.insert(std::pair<TimePoint,std::shared_ptr<const MeasurementBase>>(t,meas));
    if(!ret.second){
      std::cout << "Error: measurement already exists!" << std::endl;
    }
  }
}

void MeasurementTimeline::removeProcessedFirst(){
  assert(measMap_.size() > 0);
  lastProcessedTime_ = measMap_.begin()->first;
  measMap_.erase(measMap_.begin());
  hasProcessedTime_ = true;
}

void MeasurementTimeline::removeProcessedMeas(const TimePoint& t){
  assert(measMap_.count(t) > 0);
  measMap_.erase(t);
  lastProcessedTime_ = t;
  hasProcessedTime_ = true;
}

void MeasurementTimeline::clear(){
  measMap_.clear();
  hasProcessedTime_ = false;
}

bool MeasurementTimeline::getLastTime(TimePoint& lastTime) const{
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

TimePoint MeasurementTimeline::getMaximalUpdateTime(const TimePoint& currentTime) const{
  TimePoint maximalUpdateTime = currentTime-maxWaitTime_;
  if(!measMap_.empty()){
    maximalUpdateTime = std::max(maximalUpdateTime,measMap_.rbegin()->first+minWaitTime_);
  } else if(hasProcessedTime_){
    maximalUpdateTime = std::max(maximalUpdateTime,lastProcessedTime_+minWaitTime_);
  }
  return maximalUpdateTime;
}

void MeasurementTimeline::addAllInRange(std::set<TimePoint>& times, const TimePoint& start, const TimePoint& end) const{
  auto it = measMap_.upper_bound(start);
  while (it != measMap_.end() && it->first<= end){
    times.insert(it->first);
    ++it;
  }
}

void MeasurementTimeline::addLastInRange(std::set<TimePoint>& times, const TimePoint& start, const TimePoint& end) const{
  auto it = measMap_.upper_bound(end);
  if(it!=measMap_.begin()){
    --it;
    if(it->first > start){
      times.insert(it->first);
    }
  }
}

void MeasurementTimeline::splitMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, std::shared_ptr<const BinaryResidualBase>& res){
  addMeas(std::shared_ptr<const MeasurementBase>(),t1);
  res->splitMeasurements(measMap_.at(t2),t0,t1,t2,measMap_.at(t1),measMap_.at(t2));
}

void MeasurementTimeline::mergeMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, std::shared_ptr<const BinaryResidualBase>& res){
  res->mergeMeasurements(measMap_.at(t1),measMap_.at(t2),t0,t1,t2,measMap_.at(t2));
}

}
