/*
 * Measurement.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_
#include <map>
#include <set>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

class BinaryResidualBase;

class MeasurementBase: public State{
 public:
  MeasurementBase(const std::shared_ptr<const StateDefinition>& def): State(def){
  }
  virtual ~MeasurementBase(){};
};

class MeasurementTimeline{
 public:
  MeasurementTimeline(const bool ignoreFirst, const Duration& maxWaitTime = fromSec(0.1), const Duration& minWaitTime = Duration::zero());
  virtual ~MeasurementTimeline();
  void addMeas(const std::shared_ptr<const MeasurementBase>& meas, const TimePoint& t);
  void removeProcessedFirst();
  void removeProcessedMeas(const TimePoint& t);
  void reset();
  TimePoint getLastTime() const;
  TimePoint getMaximalUpdateTime(const TimePoint& currentTime) const;
  void addAllInRange(std::set<TimePoint>& times, const TimePoint& start, const TimePoint& end) const;
  void addLastInRange(std::set<TimePoint>& times, const TimePoint& start, const TimePoint& end) const;
  void split(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, const std::shared_ptr<const BinaryResidualBase>& res);
  void split(const std::set<TimePoint>& times, const std::shared_ptr<const BinaryResidualBase>& res);
  void merge(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, const std::shared_ptr<const BinaryResidualBase>& res);
  void mergeUndesired(const std::set<TimePoint>& times, const std::shared_ptr<const BinaryResidualBase>& res);
  void print(const TimePoint& start = TimePoint::min()) const;
 protected:
  std::map<TimePoint,std::shared_ptr<const MeasurementBase>> measMap_;
  Duration maxWaitTime_;
  Duration minWaitTime_;
  TimePoint lastProcessedTime_;
  bool ignoreFirst_;
};

}

#endif /* GIF_MEASUREMENT_HPP_ */
