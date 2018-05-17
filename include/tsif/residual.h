#ifndef TSIF_RESIDUAL_H_
#define TSIF_RESIDUAL_H_

#include "tsif/utils/common.h"
#include "tsif/model.h"

namespace tsif{

template<typename Out, typename Pre, typename Cur, typename Meas>
class Residual: public Model<Residual<Out,Pre,Cur,Meas>,Out,Pre,Cur>{
 public:
  typedef Out Output;
  typedef Pre Previous;
  typedef Cur Current;
  typedef Meas Measurement;
  std::shared_ptr<const Meas> meas_;
  double dt_;
  double w_;
  const bool isSplitable_;  // Can measurements be split into two
  const bool isMergeable_;  // Can two measurements be merged into one (should be same as isSplitable)
  const bool isMandatory_;   // Is this measurement required at every timestep (should then typically be splitable)
  bool isActive_;           // Temporary, is a measurement currently available
  Residual(bool isSplitable = true,bool isMergeable = true,bool isMandatory = true):
      isSplitable_(isSplitable),
      isMergeable_(isMergeable),
      isMandatory_(isMandatory),
      meas_(nullptr){
    dt_ = 0.1;
    w_ = 1.0;
    isActive_ = false;
    std::shared_ptr<Meas> meas(new Meas());
    meas->SetRandom();
    meas_ = meas;
  }
  virtual ~Residual(){};
  virtual int EvalRes(typename Out::Ref out, const typename Pre::CRef pre, const typename Cur::CRef cur){
    return 1;
  }
  virtual int JacPre(MatRefX J, const typename Pre::CRef pre, const typename Cur::CRef cur){
    return 1;
  }
  virtual int JacCur(MatRefX J, const typename Pre::CRef pre, const typename Cur::CRef cur){
    return 1;
  }
  int EvalImpl(typename Out::Ref out, const std::tuple<typename Pre::CRef,typename Cur::CRef> ins){
    return EvalRes(out,std::get<0>(ins),std::get<1>(ins));
  }
  template<int N, typename std::enable_if<N == 0>::type* = nullptr>
  int JacImpl(MatRefX J, const std::tuple<typename Pre::CRef,typename Cur::CRef> ins){
    return JacPre(J,std::get<0>(ins),std::get<1>(ins));
  }
  template<int N, typename std::enable_if<N == 1>::type* = nullptr>
  int JacImpl(MatRefX J, const std::tuple<typename Pre::CRef,typename Cur::CRef> ins){
    return JacCur(J,std::get<0>(ins),std::get<1>(ins));
  }
  int JacPreTest(double th, double d, const typename Pre::CRef pre, const typename Cur::CRef cur){
    return this->template JacTest<0>(th,d,std::forward_as_tuple(pre,cur));
  }
  int JacCurTest(double th, double d, const typename Pre::CRef pre, const typename Cur::CRef cur){
    return this->template JacTest<1>(th,d,std::forward_as_tuple(pre,cur));
  }
  int JacPreTest(double th, double d){
    Pre pre;
    pre.SetRandom();
    Cur cur;
    cur.SetRandom();
    return JacPreTest(th,d,pre,cur);
  }
  int JacCurTest(double th, double d){
    Pre pre;
    pre.SetRandom();
    Cur cur;
    cur.SetRandom();
    return JacCurTest(th,d,pre,cur);
  }
  virtual void SplitMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                                 std::shared_ptr<const Meas>& m0,
                                 std::shared_ptr<const Meas>& m1,
                                 std::shared_ptr<const Meas>& m2) const{
    if (isSplitable_){
      m1 = m2;
    } else{
      assert(false);
    }
  }
  virtual void MergeMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                                 std::shared_ptr<const Meas>& m0,
                                 std::shared_ptr<const Meas>& m1,
                                 std::shared_ptr<const Meas>& m2) const{
    if (isMergeable_){
      std::shared_ptr<Meas> newMeas(new Meas());
      Vec<Meas::Dim()> dif;
      m1->Boxminus(*m2,dif);
      m2->Boxplus(toSec(t1 - t0) / toSec(t2 - t0) * dif, *newMeas);
      m2 = newMeas;
    } else{
      assert(false);
    }
  }
  virtual void AddNoise(typename Out::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    AddWeight(GetWeight(),out,J_pre,J_cur);
  }
  void AddWeight(double w, typename Out::Ref out, MatRefX J_pre, MatRefX J_cur){
    out.Scale(w);
    J_pre *= w;
    J_cur *= w;
  }
  virtual double GetWeight(){
    return w_;
  }
  template<int OUT,int STA, typename std::enable_if<(STA>=0 & OUT>=0)>::type* = nullptr>
  void SetJacCur(MatRefX J, const typename Current::CRef cur, MatCRef<Output::template GetElementDim<OUT>(),Current::template GetElementDim<STA>()> Jsub){
    J.block<Output::template GetElementDim<OUT>(),Current::template GetElementDim<STA>()>(
        Output::Start(OUT),cur.Start(STA)) = Jsub;
  }
  template<int OUT,int STA, typename std::enable_if<(STA<0 | OUT<0)>::type* = nullptr>
  void SetJacCur(MatRefX J, const typename Current::CRef cur, MatCRefX Jsub){
  }
  template<int OUT,int STA, typename std::enable_if<(STA>=0 & OUT>=0)>::type* = nullptr>
  void SetJacPre(MatRefX J, const typename Previous::CRef pre, MatCRef<Output::template GetElementDim<OUT>(),Previous::template GetElementDim<STA>()> Jsub){
    J.block<Output::template GetElementDim<OUT>(),Previous::template GetElementDim<STA>()>(
        Output::Start(OUT),pre.Start(STA)) = Jsub;
  }
  template<int OUT,int STA, typename std::enable_if<(STA<0 | OUT<0)>::type* = nullptr>
  void SetJacPre(MatRefX J, const typename Previous::CRef pre, MatCRefX Jsub){
  }
};

} // namespace tsif

#endif  // TSIF_RESIDUAL_H_
