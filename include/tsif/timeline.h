#ifndef TSIF_TIMELINE_H_
#define TSIF_TIMELINE_H_

#include "element_vector.h"
#include "tsif/utils/common.h"

namespace tsif{

template<typename Measurement>
class Timeline{
 private:
  std::map<TimePoint,std::shared_ptr<const Measurement>> mm_;
  void RemoveFirst(){
    if(mm_.size() > 0){
      mm_.erase(mm_.begin());
    }
  }
 public:
  Timeline(Duration max_wait_time,Duration min_wait_time):
      max_wait_time_(max_wait_time),min_wait_time_(min_wait_time){
  }
  void Add(TimePoint t,std::shared_ptr<const Measurement> m){
    TSIF_LOG("Add measurement at " << tsif::Print(t));
    TSIF_LOGWIF(mm_.count(t) > 0, "Entry already exists for measurement!");
    mm_[t] = m;
  }
  bool HasMeas(TimePoint t){
    return mm_.count(t) > 0;
  }
  std::shared_ptr<const Measurement> Get(TimePoint t){
    return mm_.at(t);
  }
  void Clean(TimePoint t){
    while (CountSmallerOrEqual(t) > 1){ // Leave at least one measurement
      RemoveFirst();
    }
  }
  void Clear(){
    mm_.clear();
  }
  int CountSmallerOrEqual(TimePoint t){
    int count = 0;
    auto it = mm_.upper_bound(t);
    while(it != mm_.begin()){
      count++;
      it--;
    }
    return count;
  }
  TimePoint GetLastTime() const{
    if(mm_.empty()){
      return TimePoint::min();
    } else{
      return mm_.rbegin()->first;
    }
  }
  TimePoint GetFirstTime() const{
    if(mm_.empty()){
      return TimePoint::max();
    } else{
      return mm_.begin()->first;
    }
  }
  TimePoint GetMaximalUpdateTime(TimePoint current) const{
    return std::max(current-max_wait_time_,GetLastTime()+min_wait_time_);
  }
  void GetAllInRange(std::set<TimePoint>& times, TimePoint start, TimePoint end) const{
    auto it = mm_.upper_bound(start);
    while (it != mm_.end() && it->first <= end){
      times.insert(it->first);
      ++it;
    }
  }
  template<typename Residual>
  void Split(TimePoint t0, TimePoint t1, TimePoint t2, const Residual& res){
    assert(t0 <= t1 && t1 <= t2);
    TSIF_LOG("Split measurement at " << tsif::Print(t1));
    Add(t1,std::make_shared<Measurement>());
    res.SplitMeasurements(t0, t1, t2, mm_.at(t0), mm_.at(t1), mm_.at(t2));
  }
  template<typename Residual>
  void Merge(TimePoint t0, TimePoint t1, TimePoint t2, const Residual& res){
    assert(t0 <= t1 && t1 <= t2);
    TSIF_LOG("Merging measurement at " << tsif::Print(t1));
    res.MergeMeasurements(t0, t1, t2, mm_.at(t0), mm_.at(t1), mm_.at(t2));
    mm_.erase(t1);
  }
  template<typename Residual>
  void SplitAndMerge(TimePoint t0, const std::set<TimePoint>& times,const Residual& res){
    if (times.size() == 0){
      return;
    }
    if (res.isSplitable_){
      for (const auto& t : times){
        auto it = mm_.lower_bound(t);
        if (it == mm_.end()){
          TSIF_LOGW("Partial splitting only!");
          break;
        }
        if (it->first == t){
          continue; // Measurement already available
        }
        assert(it != mm_.begin());
        TimePoint previous = std::prev(it)->first;
        Split(previous, t, it->first, res);
      }
    }
    if (res.isMergeable_){
      for (auto it = mm_.begin(); it != mm_.end();){
        if (it->first <= t0 || times.count(it->first) > 0){
          ++it;
          continue; // Ignore the first or if in times
        }
        if (it->first > *times.rbegin()){
          break; // Ignore the last
        }
        if (std::next(it) == mm_.end()){
          TSIF_LOGW("Partial merging only!");
          break;
        }
        assert(it != mm_.begin());
        TimePoint previous = std::prev(it)->first;
        ++it;  // Needs to be increment before erase
        Merge(previous, std::prev(it)->first, it->first, res);
      }
    }
  }

  std::string Print(const TimePoint& start, int start_offset, double resolution) const{
    std::ostringstream out;
    const int width = mm_.empty() ? start_offset : start_offset
                      + std::max(0,(int)(std::ceil(toSec(mm_.rbegin()->first - start) / resolution)) + 1);
    std::vector<int> counts(width, 0);
    for (auto it = mm_.begin(); it != mm_.end(); ++it){
      const int x = start_offset + ceil(toSec(it->first - start) / resolution);
      if (x >= 0){
        counts.at(x)++;
      }
    }
    for (const auto& c : counts){
      if (c == 0){
        out << "-";
      } else{
        out << c;
      }
    };
    return out.str();
  }

  void SetMaxWaitTime(const double& max_wait_time){
    max_wait_time_ = fromSec(max_wait_time);
  }
  void SetMinWaitTime(const double& min_wait_time){
    min_wait_time_ = fromSec(min_wait_time);
  }

 protected:
  Duration max_wait_time_;
  Duration min_wait_time_;
};

template<>
class Timeline<MeasEmpty>{
 public:
  Timeline(Duration max_wait_time,Duration min_wait_time):
      max_wait_time_(max_wait_time),min_wait_time_(min_wait_time){
  }
  void Add(TimePoint t,std::shared_ptr<const MeasEmpty> m){
    TSIF_LOGW("Unnecessary addition of empty measurement!");
  }
  bool HasMeas(TimePoint t){
    return true;
  }
  std::shared_ptr<const MeasEmpty> Get(TimePoint t){
    return std::make_shared<MeasEmpty>();
  }
  void Clean(TimePoint t){
  }
  void Clear(){
  }
  TimePoint GetLastTime() const{
    return TimePoint::max();
  }
  TimePoint GetFirstTime() const{
    return TimePoint::min();
  }
  TimePoint GetMaximalUpdateTime(TimePoint current) const{
    return TimePoint::max();
  }
  void GetAllInRange(std::set<TimePoint>& times, TimePoint start, TimePoint end) const{
  }
  template<typename Residual>
  void SplitAndMerge(TimePoint time, const std::set<TimePoint>& times,const Residual& res){
  }
  std::string Print(const TimePoint& start, int start_offset, double resolution) const{
    const int c = start_offset > 4 ? start_offset - 4 : 0;
    std::ostringstream out;
    for(int i=0;i<c;i++) out << "-";
    out << "MeasEmpty";
    for(int i=0;i<c;i++) out << "-";
    return out.str();
  }

  void SetMaxWaitTime(const double& max_wait_time){
    max_wait_time_ = fromSec(max_wait_time);
  }
  void SetMinWaitTime(const double& min_wait_time){
    min_wait_time_ = fromSec(min_wait_time);
  }

 protected:
  Duration max_wait_time_;
  Duration min_wait_time_;
};

} // namespace tsif

#endif  // TSIF_TIMELINE_H_
