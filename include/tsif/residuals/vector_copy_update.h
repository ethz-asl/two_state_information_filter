#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int N>
class VecMeas: public ElementVector<Element<Vec<N>,0>>{
 public:
  VecMeas(): ElementVector<Element<Vec<N>,0>>(Vec<N>::Zero()){}
  VecMeas(const Vec<N>& vec): ElementVector<Element<Vec<N>,0>>(vec){}
  const Vec<N>& GetMeas() const{
    return this->template Get<0>();
  }
  Vec<N>& GetMeas(){
    return this->template Get<0>();
  }
};

template<int Y, int STA, int N>
using VectorCopyUpdateBase = Residual<ElementVector<Element<Vec<N>,Y>>,
                                                    ElementVector<>,
                                                    ElementVector<Element<Vec<N>,STA>>,
                                                    VecMeas<N>>;

template<int STA, int N, bool SNM>
class VectorCopyUpdate: public VectorCopyUpdateBase<0,STA,N>{
 public:
  typedef VectorCopyUpdateBase<0,STA,N> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  double huber_threshold_;
  VectorCopyUpdate(): Base(SNM,SNM,true), huber_threshold_(1000.){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = cur.template Get<STA>() - meas_->GetMeas();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    this->template SetJacCur<0,STA>(J,cur,Mat<N>::Identity());
    return 0;
  }
  virtual void AddNoise(typename Output::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    //compute weight using huber loss function
    double weight = w_;
    const double norm = out.template Get<0>().norm();
    if(norm > huber_threshold_) weight *= sqrt(huber_threshold_ * (norm - 0.5 * huber_threshold_))/norm;
    //scale the innovation and jacobians
    out.Scale(weight);
    J_pre *= weight;
    J_cur *= weight;
  }
};

} // namespace tsif
