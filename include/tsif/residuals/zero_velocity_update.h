#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class ZeroVelocityTrigger: public ElementVector<Element<bool,0>>{
 public:
  ZeroVelocityTrigger(): ElementVector<Element<bool,0>>(false){}
  ZeroVelocityTrigger(const bool& trigger): ElementVector<Element<bool,0>>(trigger){}
  const bool& GetTrigger() const{
    return Get<0>();
  }
  bool& GetTrigger(){
    return Get<0>();
  }
};

template<int Y_V, int Y_OMEGA, int B_V_IB, int B_OMEGA_IB>
using ZeroVelocityUpdateBase = Residual<ElementVector<Element<Vec3,Y_V>, Element<Vec3,Y_OMEGA>>,
                                        ElementVector<>,
                                        ElementVector<Element<Vec3,B_V_IB>,Element<Vec3,B_OMEGA_IB>>,
                                        ZeroVelocityTrigger>;

template<int B_V_IB, int B_OMEGA_IB>
class ZeroVelocityUpdate: public ZeroVelocityUpdateBase<0,1,B_V_IB, B_OMEGA_IB>{
 public:
  typedef ZeroVelocityUpdateBase<0,1,B_V_IB, B_OMEGA_IB> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  ZeroVelocityUpdate(): Base(false,false,false) {}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    if (meas_->GetTrigger()){
      out.template Get<0>() = cur.template Get<B_V_IB>();
      out.template Get<1>() = cur.template Get<B_OMEGA_IB>();
    }
    else {
      out.template Get<0>() = Vec3::Zero();
      out.template Get<1>() = Vec3::Zero();
    }
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    if (meas_->GetTrigger()){
      this->template SetJacCur<0,B_V_IB>(J,cur,Mat3::Identity());
      this->template SetJacCur<1,B_OMEGA_IB>(J,cur,Mat3::Identity());
    }
    else {
      this->template SetJacCur<0,B_V_IB>(J,cur,Mat3::Zero());
      this->template SetJacCur<1,B_OMEGA_IB>(J,cur,Mat3::Zero());
    }
    return 0;
  }
};

} // namespace tsif
