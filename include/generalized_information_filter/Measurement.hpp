#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_
#include <map>
#include <set>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF {

class BinaryResidualBase;

class MeasurementTimeline {
 public:
  MeasurementTimeline(const bool ignoreFirst, const Duration& maxWaitTime =
                          fromSec(0.1),
                      const Duration& minWaitTime = Duration::zero());
  virtual ~MeasurementTimeline();
  void addMeas(const std::shared_ptr<const StateBase>& meas,
               const TimePoint& t);
  bool getMeas(const TimePoint& t, std::shared_ptr<const StateBase>& meas);
  void removeProcessedFirst();
  void removeProcessedMeas(const TimePoint& t);
  void reset();
  TimePoint getLastTime() const;
  TimePoint getMaximalUpdateTime(const TimePoint& currentTime) const;
  void getAllInRange(std::set<TimePoint>& times, const TimePoint& start,
                     const TimePoint& end) const;
  void getLastInRange(std::set<TimePoint>& times, const TimePoint& start,
                      const TimePoint& end) const;
  void split(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const std::shared_ptr<const BinaryResidualBase>& res);
  void split(const std::set<TimePoint>& times,
             const std::shared_ptr<const BinaryResidualBase>& res);
  void merge(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const std::shared_ptr<const BinaryResidualBase>& res);
  void mergeUndesired(const std::set<TimePoint>& times,
                      const std::shared_ptr<const BinaryResidualBase>& res);
  void removeOutdated(const TimePoint& time);
  void print(const TimePoint& start = TimePoint::min(), int startOffset = 0,
             double resolution = 0.01) const;
 protected:
  std::map<TimePoint, std::shared_ptr<const StateBase>> measMap_;
  Duration maxWaitTime_;
  Duration minWaitTime_;
  TimePoint lastProcessedTime_;
  bool ignoreFirst_;
};

}
#endif /* GIF_MEASUREMENT_HPP_ */
