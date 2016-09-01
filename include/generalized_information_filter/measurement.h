#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_
#include <map>
#include <set>

#include "element-vector.h"
#include "generalized_information_filter/common.h"

namespace GIF {

class BinaryResidualBase;

class MeasurementTimeline {
 public:
  MeasurementTimeline(const bool ignoreFirst, const Duration& maxWaitTime =
                          fromSec(0.1),
                      const Duration& minWaitTime = Duration::zero());
  virtual ~MeasurementTimeline();
  void AddMeasurement(const std::shared_ptr<const ElementVectorBase>& meas,
               const TimePoint& t);
  bool GetMeasurement(const TimePoint& t, std::shared_ptr<const ElementVectorBase>& meas);
  void RemoveProcessedFirst();
  void RemoveProcessedMeas(const TimePoint& t);
  void Reset();
  TimePoint GetLastTime() const;
  TimePoint GetMaximalUpdateTime(const TimePoint& currentTime) const;
  void GetAllInRange(std::set<TimePoint>& times, const TimePoint& start,
                     const TimePoint& end) const;
  void GetLastInRange(std::set<TimePoint>& times, const TimePoint& start,
                      const TimePoint& end) const;
  void Split(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const BinaryResidualBase* res);
  void Split(const std::set<TimePoint>& times,
             const BinaryResidualBase* res);
  void Merge(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const BinaryResidualBase* res);
  void MergeUndesired(const std::set<TimePoint>& times,
                      const BinaryResidualBase* res);
  void RemoveOutdated(const TimePoint& time);
  void Print(const TimePoint& start = TimePoint::min(), int startOffset = 0,
             double resolution = 0.01) const;
 protected:
  std::map<TimePoint, std::shared_ptr<const ElementVectorBase>> measMap_;
  Duration maxWaitTime_;
  Duration minWaitTime_;
  TimePoint lastProcessedTime_;
  bool ignoreFirst_;
};

}
#endif /* GIF_MEASUREMENT_HPP_ */
