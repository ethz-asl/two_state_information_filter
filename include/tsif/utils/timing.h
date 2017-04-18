#ifndef TSIF_TIMING_HPP_
#define TSIF_TIMING_HPP_

#include <iostream>
#include <chrono>

namespace tsif{

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double sec){
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration){
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}
static std::string Print(TimePoint t){
  std::ostringstream out;
  out.precision(15);
  out << ((double)t.time_since_epoch().count()*Duration::period::num)/(Duration::period::den);
  return out.str();
}

class Timer{
 public:
  Timer(){
    start_ = Clock::now();
    last_ = start_;
  }
  double GetIncr(){
    TimePoint now = Clock::now();
    double incr = toSec(now-last_);
    last_ = now;
    return incr;
  }
  double GetFull(){
    last_ = Clock::now();
    return toSec(last_-start_);
  }
  TimePoint start_;
  TimePoint last_;
};

} // namespace tsif

#endif /* TSIF_TIMING_HPP_ */
