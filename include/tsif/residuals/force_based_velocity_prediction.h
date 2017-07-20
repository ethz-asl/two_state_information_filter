#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{
class MeasContactForces: public ElementVector<Element<Vec3,0>>{
 public:
  MeasContactForces(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasContactForces(const Vec3& forces): ElementVector<Element<Vec3,0>>(forces){}
  const Vec3& GetMeas() const{
    return Get<0>();
  }
  Vec3& GetMeas(){
    return Get<0>();
  }
};

template<int Y, int PHI_IB, int B_V_IB, int B_OMEGA_IB, int M_R>
using ForceBasedVelocityPredictionBase = Residual<ElementVector<Element<Vec3,Y>>,
                                                  ElementVector<Element<Quat,PHI_IB>,Element<Vec3,B_V_IB>,Element<Vec3,B_OMEGA_IB>,Element<double,M_R>>,
                                                  ElementVector<Element<Vec3,B_V_IB>>,
                                                  MeasContactForces>;

template<int PHI_IB, int B_V_IB, int B_OMEGA_IB, int M_R>
class ForceBasedVelocityPrediction: public ForceBasedVelocityPredictionBase<0, PHI_IB, B_V_IB, B_OMEGA_IB, M_R>{
 public:
  typedef ForceBasedVelocityPredictionBase<0, PHI_IB, B_V_IB, B_OMEGA_IB, M_R> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  Vec3 g_;
  ForceBasedVelocityPrediction(): Base(true,true,false), g_(0,0,-9.81){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = cur.template Get<B_V_IB>()
        - (Mat3::Identity() - SSM(dt_*pre.template Get<B_OMEGA_IB>()))*pre.template Get<B_V_IB>()
        - dt_*(meas_->GetMeas()/pre.template Get<M_R>()+pre.template Get<PHI_IB>().inverse().toRotationMatrix()*g_);
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    this->template SetJacPre<0,PHI_IB>(J,pre,-pre.template Get<PHI_IB>().inverse().toRotationMatrix()*SSM(g_*dt_));
    this->template SetJacPre<0,B_V_IB>(J,pre, -(Mat3::Identity()-SSM(dt_*pre.template Get<B_OMEGA_IB>())));
    this->template SetJacPre<0,B_OMEGA_IB>(J,pre,-SSM(dt_*pre.template Get<B_V_IB>()));
    this->template SetJacPre<0,M_R>(J,pre,dt_*meas_->GetMeas()/(pre.template Get<M_R>()*pre.template Get<M_R>()));
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    this->template SetJacCur<0,B_V_IB>(J,cur,Mat3::Identity());
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif
