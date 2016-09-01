#include "generalized_information_filter/measurement.h"
#include "generalized_information_filter/binary-residual.h"

namespace GIF {

MeasurementTimeline::MeasurementTimeline(const bool ignoreFirst,
                                         const Duration& maxWaitTime,
                                         const Duration& minWaitTime) {
  ignoreFirst_ = ignoreFirst;
  maxWaitTime_ = maxWaitTime;
  minWaitTime_ = minWaitTime;
  lastProcessedTime_ = TimePoint::min();
}

MeasurementTimeline::~MeasurementTimeline() {
}

void MeasurementTimeline::AddMeasurement(const std::shared_ptr<const ElementVectorBase>& meas,
                                  const TimePoint& t) {
  // Discard first measurement in binary case
  if (ignoreFirst_ && lastProcessedTime_ == TimePoint::min()) {
    lastProcessedTime_ = t;
    return;
  }
  if (t <= lastProcessedTime_) {
    std::cout
        << "Error: adding measurements before last processed time (will be "
        << "discarded)" << std::endl;
  } else {
    std::pair<std::map<TimePoint, std::shared_ptr<const ElementVectorBase>>::iterator,
        bool> ret;
    ret = measMap_.insert(
        std::pair<TimePoint, std::shared_ptr<const ElementVectorBase>>(t, meas));
    if (!ret.second) {
      std::cout << "Error: measurement already exists!" << std::endl;
    }
  }
}

bool MeasurementTimeline::GetMeasurement(const TimePoint& t,
                                  std::shared_ptr<const ElementVectorBase>& meas) {
  auto it = measMap_.find(t);
  if (it == measMap_.end()) {
    return false;
  } else {
    meas = it->second;
    return true;
  }
}

void MeasurementTimeline::RemoveProcessedFirst() {
  assert(measMap_.size() > 0);
  lastProcessedTime_ = measMap_.begin()->first;
  measMap_.erase(measMap_.begin());
}

void MeasurementTimeline::RemoveProcessedMeas(const TimePoint& t) {
  assert(measMap_.count(t) > 0);
  measMap_.erase(t);
  lastProcessedTime_ = t;
}

void MeasurementTimeline::Reset() {
  measMap_.clear();
  lastProcessedTime_ = TimePoint::min();
}

TimePoint MeasurementTimeline::GetLastTime() const {
  if (!measMap_.empty()) {
    return measMap_.rbegin()->first;
  } else {
    return lastProcessedTime_;
  }
}

TimePoint MeasurementTimeline::GetMaximalUpdateTime(
    const TimePoint& currentTime) const {
  TimePoint maximalUpdateTime = currentTime - maxWaitTime_;
  if (!measMap_.empty()) {
    maximalUpdateTime = std::max(maximalUpdateTime,
                                 measMap_.rbegin()->first + minWaitTime_);
  } else {
    maximalUpdateTime = std::max(maximalUpdateTime,
                                 lastProcessedTime_ + minWaitTime_);
  }
  return maximalUpdateTime;
}

void MeasurementTimeline::GetAllInRange(std::set<TimePoint>& times,
                                        const TimePoint& start,
                                        const TimePoint& end) const {
  auto it = measMap_.upper_bound(start);
  while (it != measMap_.end() && it->first <= end) {
    times.insert(it->first);
    ++it;
  }
}

void MeasurementTimeline::GetLastInRange(std::set<TimePoint>& times,
                                         const TimePoint& start,
                                         const TimePoint& end) const {
  auto it = measMap_.upper_bound(end);
  if (it != measMap_.begin()) {
    --it;
    if (it->first > start) {
      times.insert(it->first);
    }
  }
}

void MeasurementTimeline::Split(
    const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
    const BinaryResidualBase* res) {
  assert(t0 <= t1 && t1 <= t2);
  AddMeasurement(std::shared_ptr<const ElementVectorBase>(), t1);
  res->SplitMeasurements(t0, t1, t2, measMap_.at(t2), measMap_.at(t1),
                         measMap_.at(t2));
}

void MeasurementTimeline::Split(
    const std::set<TimePoint>& times,
    const BinaryResidualBase* res) {
  for (const auto& t : times) {
    auto it = measMap_.lower_bound(t);
    if (it == measMap_.end()) {
      std::cout << "Error: range error while splitting!" << std::endl;
      continue;
    }
    if (it->first == t) {
      // Measurement already available
      continue;
    }
    TimePoint previous =
        (it == measMap_.begin()) ? lastProcessedTime_ : std::prev(it)->first;
    Split(previous, t, it->first, res);
  }
}

void MeasurementTimeline::Merge(
    const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
    const BinaryResidualBase* res) {
  assert(t0 <= t1 && t1 <= t2);
  res->MergeMeasurements(t0, t1, t2, measMap_.at(t1), measMap_.at(t2), measMap_.at(t2));
  measMap_.erase(t1);  // does not count as processed
}

void MeasurementTimeline::MergeUndesired(
    const std::set<TimePoint>& times,
    const BinaryResidualBase* res) {
  // Merge measurements such that only timepoints remain which are in times or
  // past its end.
  if (times.size() == 0) {
    return;
  }
  for (auto it = measMap_.begin(); it != measMap_.end();) {
    if (it->first > *times.rbegin()) {
      break;
    }
    if (times.count(it->first) > 0) {
      ++it;
      continue;
    }
    if (next(it) == measMap_.end()) {
      std::cout << "Error: range error while merging!" << std::endl;
      break;
    }
    TimePoint previous =
        (it == measMap_.begin()) ? lastProcessedTime_ : std::prev(it)->first;
    ++it;  // Needs to be increment before erase
    Merge(previous, std::prev(it)->first, it->first, res);
  }
}

void MeasurementTimeline::RemoveOutdated(const TimePoint& time) {
  while (!measMap_.empty() && measMap_.begin()->first <= time) {
    RemoveProcessedFirst();
  }
}

void MeasurementTimeline::Print(const TimePoint& start, int startOffset,
                                double resolution) const {
  const int width =
      measMap_.empty() ?
          startOffset :
          startOffset
              + ceil(toSec(measMap_.rbegin()->first - start) / resolution) + 1;
  std::vector<int> counts(width, 0);
  for (auto it = measMap_.begin(); it != measMap_.end(); ++it) {
    const int x = startOffset + ceil(toSec(it->first - start) / resolution);
    if (x >= 0) {
      counts[x]++;
    }
  }
  for (const auto& c : counts) {
    if (c == 0) {
      std::cout << "-";
    } else {
      std::cout << c;
    }
  }
  std::cout << std::endl;
}

}
