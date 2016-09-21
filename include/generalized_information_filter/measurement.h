#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_
#include <map>
#include <set>

#include "element-vector.h"
#include "generalized_information_filter/common.h"

namespace GIF {

class EmptyMeas : public ElementVector {
 public:
  EmptyMeas(): ElementVector(std::make_shared<ElementVectorDefinition>()){};
};

class BinaryResidualBase;

/*! \brief MeasurementTimeline
 *         Class for managing measurements in a timeline. Provides various helpers.
 *         Assume: last_processed_time_ is always smaller than all other measurement times.
 */
class MeasurementTimeline {
 public:
  MeasurementTimeline(const bool ignore_first, const Duration& max_wait_time,
                      const Duration& min_wait_time);
  virtual ~MeasurementTimeline();
  void AddMeasurement(const std::shared_ptr<const ElementVectorBase>& meas, const TimePoint& t);
  bool GetMeasurement(const TimePoint& t, std::shared_ptr<const ElementVectorBase>& meas);
  void RemoveProcessedFirst();
  void Reset();
  TimePoint GetLastTime() const;
  TimePoint GetFirstTime() const;
  TimePoint GetMaximalUpdateTime(const TimePoint& current_time) const;
  void GetAllInRange(std::set<TimePoint>& times, const TimePoint& start,
                     const TimePoint& end) const;
  void GetLastInRange(std::set<TimePoint>& times, const TimePoint& start,
                      const TimePoint& end) const;
  void Split(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const BinaryResidualBase* res);
  void Split(const std::set<TimePoint>& times, const BinaryResidualBase* res);
  void Merge(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
             const BinaryResidualBase* res);
  void MergeUndesired(const std::set<TimePoint>& times, const BinaryResidualBase* res);
  void RemoveOutdated(const TimePoint& time);
  void Print(const TimePoint& start, int start_offset, double resolution) const;
 protected:
  std::map<TimePoint, std::shared_ptr<const ElementVectorBase>> meas_map_;
  Duration max_wait_time_;
  Duration min_wait_time_; // Should always be zero for binary residual
  TimePoint last_processed_time_;
  bool drop_first_;
};

}
#endif /* GIF_MEASUREMENT_HPP_ */
