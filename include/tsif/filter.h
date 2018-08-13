#ifndef TSIF_FILTER_H_
#define TSIF_FILTER_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"
#include "tsif/timeline.h"

namespace tsif{

template<typename... Residuals>
class Filter{
 public:
  static constexpr int kN = sizeof...(Residuals);
  typedef typename MergeTrait<typename Residuals::Previous...,typename Residuals::Current...>::Type State;
  typedef std::tuple<Residuals...> ResidualTuple;
  typedef std::tuple<Timeline<typename Residuals::Measurement>...> TimelineTuple;
  alignas(16) ResidualTuple residuals_;
  alignas(16) TimelineTuple timelines_;

  Filter(): timelines_(Timeline<typename Residuals::Measurement>(fromSec(0.1),fromSec(0.0))...),
            I_(State::Dim(),State::Dim()){
    is_initialized_ = false;
    include_max_ = false;
    max_iter_ = 1;
    th_iter_ = 0.1;
    iter_ = 0;
    weightedDelta_ = 0;
  }
  virtual ~Filter(){}

  template<int N>
  void AddMeas(TimePoint t,std::shared_ptr<const typename std::tuple_element<N,ResidualTuple>::type::Measurement> m){
    TSIF_LOGWIF(t < time_,"Warning: adding measurement in past!");
    TSIF_LOGWIF(is_initialized_ && (t-time_ > fromSec(10)),"Warning: measurement far in future!");
    std::get<N>(timelines_).Add(t,m);
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  TimePoint GetMaxUpdateTime(TimePoint current){
    return std::min(std::get<C>(timelines_).GetMaximalUpdateTime(current), GetMaxUpdateTime<C+1>(current));
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  TimePoint GetMaxUpdateTime(TimePoint current){
    return current;
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  TimePoint GetCurrentTime(){
    TimePoint lastTime = std::get<C>(timelines_).GetLastTime();
    return lastTime == TimePoint::max() ? GetCurrentTime<C+1>() : std::max(std::get<C>(timelines_).GetLastTime(), GetCurrentTime<C+1>());
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  TimePoint GetCurrentTime(){
    return TimePoint::min();
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  TimePoint GetMinMaxTime(){
    return std::min(std::get<C>(timelines_).GetLastTime(), GetMinMaxTime<C+1>());
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  TimePoint GetMinMaxTime(){
    return TimePoint::max();
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  TimePoint GetMaxMinTime(){
    return std::max(std::get<C>(timelines_).GetFirstTime(), GetMaxMinTime<C+1>());
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  TimePoint GetMaxMinTime(){
    return TimePoint::min();
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void PrintTimelines(const TimePoint& start, int startOffset, double resolution){
    PrintTimelines<C+1>(start,startOffset,resolution);
    TSIF_LOG(std::get<C>(timelines_).Print(start, startOffset, resolution));
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void PrintTimelines(const TimePoint& start, int startOffset, double resolution){
    std::ostringstream out;
    for (int i = 0; i < startOffset; i++){
      out << " ";
    }
    out << "|";
    TSIF_LOG(out.str());
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void GetTimeList(std::set<TimePoint>& times, TimePoint maxUpdateTime, bool includeMax) const{
    if (!std::get<C>(residuals_).isMergeable_){
      std::get<C>(timelines_).GetAllInRange(times, time_, maxUpdateTime);
    }
    GetTimeList<C+1>(times,maxUpdateTime,includeMax);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void GetTimeList(std::set<TimePoint>& times, TimePoint maxUpdateTime, bool includeMax) const{
    if (includeMax && maxUpdateTime > time_){
      times.insert(maxUpdateTime);
    }
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void SplitAndMerge(TimePoint time, const std::set<TimePoint>& times){
    std::get<C>(timelines_).SplitAndMerge(time,times,std::get<C>(residuals_));
    SplitAndMerge<C+1>(time,times);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void SplitAndMerge(TimePoint time, const std::set<TimePoint>& times){}

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void Clean(TimePoint time){
    std::get<C>(timelines_).Clean(time);
    Clean<C+1>(time);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void Clean(TimePoint time){}

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void Clear(){
    std::get<C>(timelines_).Clear();
    Clear<C+1>();
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void Clear(){}

  virtual void Init(TimePoint t){
    if(GetMinMaxTime() != TimePoint::min()){
      TSIF_LOG("Initializing state at t = " << Print(t));
      startTime_ = t;
      time_ = t;
      state_.SetIdentity();
      I_.setIdentity();
      is_initialized_ = true;
    }
  }
  virtual void ComputeLinearizationPoint(const TimePoint& t){
    curLinState_ = state_;
  }
  virtual void PreProcess(){};
  virtual void PostProcess(){};

  void Update(){
    // Initialize if possible
    if(!is_initialized_){
      Init(GetMaxMinTime());
    }

    if(is_initialized_){
      TSIF_LOG("Timelines before processing:");
      PrintTimelines(time_, 20, 0.001);
      TSIF_LOG("State time:\t" << Print(time_));
      TimePoint current = GetCurrentTime();
      TSIF_LOG("Current time:\t" << Print(current));
      TimePoint maxUpdateTime = GetMaxUpdateTime(current);
      TSIF_LOG("Maximal update time:\t" << Print(maxUpdateTime));
      std::set<TimePoint> times;
      GetTimeList(times, maxUpdateTime, include_max_);
      std::ostringstream out;
      out << "Update times:\t";
      for (const auto& t : times){
        out << Print(t) << " ";
      }
      TSIF_LOG(out.str());
      SplitAndMerge(time_,times);
      TSIF_LOG("Timelines after split and merging:");
      PrintTimelines(time_, 20, 0.001);

      // Carry out updates
      for (const auto& t : times){
        MakeUpdateStep(t);
      }
      Clean(time_);
      TSIF_LOG("Timelines after cleaning:");
      PrintTimelines(time_, 20, 0.001);
    }
  }

  void MakeUpdateStep(TimePoint t){
    // Compute linearisation point
    ComputeLinearizationPoint(t);

    // Check available measurements and prepare residuals
    int innDim = PreProcessResidual(t);
    PreProcess();

    // Temporaries
    y_.resize(innDim);
    y_.setZero();
    JacPre_.resize(innDim, State::Dim());
    JacPre_.setZero();
    JacCur_.resize(innDim, State::Dim());
    JacCur_.setZero();

    weightedDelta_ = th_iter_;
    MatX newInf(State::Dim(),State::Dim());
    for(iter_=0;iter_<max_iter_ && weightedDelta_ >= th_iter_;iter_++){
      ConstructProblem(0);
      TSIF_LOG("Innovation:\t" << y_.transpose());
      TSIF_LOG("JacPre:\n" << JacPre_);
      TSIF_LOG("JacCur:\n" << JacCur_);

      // Compute Kalman Update // TODO use more efficient form
      MatX D = I_ + JacPre_.transpose() * JacPre_;
      MatX J(innDim,innDim); J.setIdentity();
#if TSIF_VERBOSE > 0
      Eigen::JacobiSVD<MatX> svdD(D);
      const double condD = svdD.singularValues()(0) / svdD.singularValues()(svdD.singularValues().size()-1);
      TSIF_LOG("D condition number:\n" << condD);
#endif
      MatX S = JacCur_.transpose() * (J - JacPre_ * D.inverse() * JacPre_.transpose());
      newInf = S * JacCur_;
      newInf = 0.5*(newInf + newInf.transpose().eval());
      Eigen::LDLT<MatX> I_LDLT(newInf);
#if TSIF_VERBOSE > 0
      Eigen::JacobiSVD<MatX> svdI(newInf);
      const double condI = svdI.singularValues()(0) / svdI.singularValues()(svdI.singularValues().size()-1);
      TSIF_LOG("I condition number:\n" << condI);
#endif
      TSIF_LOGEIF((I_LDLT.info() != Eigen::Success),"Computation of Iinv failed");
      VecX dx = -I_LDLT.solve(S * y_);

      // Apply Kalman Update
      State newState = curLinState_;
      curLinState_.Boxplus(dx, newState);
      curLinState_ = newState;
      weightedDelta_ = sqrt((dx.dot(newInf*dx))/dx.size());
      TSIF_LOG("iter: " << iter_ << "\tw: " << sqrt((dx.dot(dx))/dx.size()) << "\twd: " << weightedDelta_);
    }
    TSIF_LOGWIF(weightedDelta_ >= th_iter_, "Reached maximal iterations:" << iter_);

    state_ = curLinState_;
    I_ = newInf;
    TSIF_LOG("State after Update:\n" << state_.Print());
    TSIF_LOG("Information matrix:\n" << I_);

    // Post Processing
    PostProcess();
    time_ = t;
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  int PreProcessResidual(TimePoint t){
    std::get<C>(residuals_).isActive_ = std::get<C>(timelines_).HasMeas(t);
    assert(std::get<C>(residuals_).isActive_ || !std::get<C>(residuals_).isMandatory_);
    if(std::get<C>(residuals_).isActive_){
      std::get<C>(residuals_).dt_ = toSec(t-time_);
      std::get<C>(residuals_).meas_ = std::get<C>(timelines_).Get(t);
    }
    return PreProcessResidual<C+1>(t) + std::tuple_element<C,ResidualTuple>::type::Output::Dim();
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  int PreProcessResidual(TimePoint t){
    return 0;
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void ConstructProblem(int start){
    typedef typename std::tuple_element<C,ResidualTuple>::type::Output Output;
    static_assert(Output::kIsVectorSpace, "Residual must be vector space!");
    if(std::get<C>(residuals_).isActive_ && Output::Dim() > 0){
      Output ySub;
      std::get<C>(residuals_).EvalRes(ySub,state_,curLinState_);
      std::get<C>(residuals_).JacPre(JacPre_.block<Output::Dim(),State::Dim()>(start,0),state_,curLinState_);
      std::get<C>(residuals_).JacCur(JacCur_.block<Output::Dim(),State::Dim()>(start,0),state_,curLinState_);
      std::get<C>(residuals_).AddNoise(ySub, JacPre_.block<Output::Dim(),State::Dim()>(start,0),
                                       JacCur_.block<Output::Dim(),State::Dim()>(start,0),
                                       state_,curLinState_);
      ySub.GetVec(y_.block<Output::Dim(),1>(start, 0));
    }
    ConstructProblem<C+1>(start+Output::Dim()*std::get<C>(residuals_).isActive_);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void ConstructProblem(int start){
  }

  template<int N = 0, typename std::enable_if<(N < kN)>::type* = nullptr>
  int JacTestAll(double th, double d, const State& pre, const State& cur){
    std::get<N>(residuals_).JacPreTest(th,d,pre,cur);
    std::get<N>(residuals_).JacCurTest(th,d,pre,cur);
    return JacTestAll<N+1>(th,d,pre,cur);
  }

  template<int N = 0, typename std::enable_if<(N >= kN)>::type* = nullptr>
  int JacTestAll(double th, double d, const State& pre, const State& cur){
    return 0;
  }

  template<int N = 0, typename std::enable_if<(N < kN)>::type* = nullptr>
  int JacTestAll(double th, double d){
    std::get<N>(residuals_).JacPreTest(th,d);
    std::get<N>(residuals_).JacCurTest(th,d);
    return JacTestAll<N+1>(th,d);
  }

  template<int N = 0, typename std::enable_if<(N >= kN)>::type* = nullptr>
  int JacTestAll(double th, double d){
    return 0;
  }

  const State& GetState() const{
    return state_;
  }
  State& GetState(){
    return state_;
  }

  template<int N = 0, int C = 0, bool B = false, typename std::enable_if<(N < kN)>::type* = nullptr>
  std::string PrintConnectivityRes(){
    std::ostringstream out;
    if(!B){
      if(std::tuple_element<N,ResidualTuple>::type::Previous::template HasId<State::template GetId<C>()>()){
        out << "--" << " ";
      }else{
        out << "   ";
      }
      if(C==State::kN-1){
        if(N<10) out << " ";
        out << N << " ";
      }
    } else {
      if(std::tuple_element<N,ResidualTuple>::type::Current::template HasId<State::template GetId<C>()>()){
        out << "--" << " ";
      }else{
        out << "   ";
      }
      if(C==State::kN-1) out << std::endl;
    }
    out << PrintConnectivityRes<N+(C==State::kN-1 & B),(C==State::kN-1)?0:C+1,((C==State::kN-1)?!B:B)>();
    return out.str();
  }
  template<int N = 0, int C = 0, bool B = false, typename std::enable_if<(N >= kN)>::type* = nullptr>
  std::string PrintConnectivityRes(){
    return "";
  }

  std::string PrintConnectivity(){
    std::ostringstream out;
    for(int i=0;i<State::kN;i++){
      if(i<10) out << " ";
      out << i << " ";
    }
    out << "   ";
    for(int i=0;i<State::kN;i++){
      if(i<10) out << " ";
      out << i << " ";
    }
    out << std::endl;
    out << PrintConnectivityRes();
    return out.str();
  }

  MatX GetCovariance() const {
    const int state_dim = State::Dim();
    MatX identity(state_dim, state_dim);
    identity.setIdentity();
    const MatX cov = I_.llt().solve(identity);
    return cov;
  }

  MatX GetInformation() const { return I_; }

  void Uninitialize() { is_initialized_ = false; }

  bool IsInitialized() const { return is_initialized_; }

  template <int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void SetMaxWaitTimes(double max_wait_time) {
    std::get<C>(timelines_).SetMaxWaitTime(max_wait_time);
    SetMaxWaitTimes<C + 1>(max_wait_time);
  }
  template <int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void SetMaxWaitTimes(double max_wait_time) {}
  
  template <int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void SetMinWaitTimes(double min_wait_time) {
    std::get<C>(timelines_).SetMinWaitTime(min_wait_time);
    SetMinWaitTimes<C + 1>(min_wait_time);
  }
  template <int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void SetMinWaitTimes(double min_wait_time) {}

 protected:
  bool is_initialized_;
  bool include_max_;
  TimePoint startTime_;
  TimePoint time_;
  State state_;
  State curLinState_;
  MatX I_;
  VecX y_;
  MatX JacPre_;
  MatX JacCur_;
  int max_iter_;
  int iter_;
  double weightedDelta_;
  double th_iter_;
};

} // namespace tsif

#endif  // TSIF_FILTER_H_
