# pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int Y, int STA>
using ScalarRandomWalkBase = Residual<ElementVector<Element<Vec<2>,Y>>,
                                      ElementVector<Element<double,STA>>,
                                      ElementVector<Element<double,STA>>,
                                      MeasEmpty>;

template<int STA>
class ScalarRandomWalk: public ScalarRandomWalkBase<0, STA>{
 public:
  typedef ScalarRandomWalkBase<0, STA> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  ScalarRandomWalk(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out. template Get<0>()(0) = pre.template Get<STA>()-cur.template Get<STA>();
    out. template Get<0>()(1) = 0.;
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    this->template SetJacCur<0,STA>(J,pre,Vec<2>(1.,0.));
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    this->template SetJacCur<0,STA>(J,cur,Vec<2>(-1.,0.));
    return 0;
  }
  virtual double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif
