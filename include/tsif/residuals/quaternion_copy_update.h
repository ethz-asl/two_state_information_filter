#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class QuatMeas: public ElementVector<Element<Quat,0>>{
 public:
  QuatMeas(): ElementVector<Element<Quat,0>>(Quat::Identity()){}
  QuatMeas(const Quat& quat): ElementVector<Element<Quat,0>>(quat){}
  const Quat& GetMeas() const{
    return Get<0>();
  }
  Quat& GetMeas(){
    return Get<0>();
  }

};

template<int Y, int STA>
using QuaternionCopyUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                                    ElementVector<>,
                                                    ElementVector<Element<Quat,STA>>,
                                                    QuatMeas>;

template<int STA, bool SNM>
class QuaternionCopyUpdate: public QuaternionCopyUpdateBase<0,STA>{
 public:
  typedef QuaternionCopyUpdateBase<0,STA> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  double huber_threshold_;
  QuaternionCopyUpdate(): Base(SNM,SNM,true), huber_threshold_(1000.){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = Log(cur.template Get<STA>()*meas_->GetMeas().inverse());
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    const Vec3 y = Log(cur.template Get<STA>()*meas_->GetMeas().inverse());
    const Mat3 C =  cur.template Get<STA>().toRotationMatrix();
    const Mat3 G_inv = GammaMat(y).inverse();
    this->template SetJacCur<0,STA>(J,cur,G_inv);
    return 0;
  }
  virtual void AddNoise(typename Output::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    //compute weight using huber loss function
    double weight = w_;
    const double norm = out.template Get<0>().norm();
    if(norm > huber_threshold_) weight *= sqrt(2*huber_threshold_ * (norm - 0.5 * huber_threshold_))/norm;
    //scale the innovation and jacobians
    out.Scale(weight);
    J_pre *= weight;
    J_cur *= weight;
  }
};

} // namespace tsif
