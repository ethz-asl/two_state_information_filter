#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasLvel: public ElementVector<Element<Vec3,0>>{
 public:
  MeasLvel(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasLvel(const Vec3& vel): ElementVector<Element<Vec3,0>>(vel){}
  const Vec3& GetMeas() const{
    return Get<0>();
  }
  Vec3& GetMeas(){
    return Get<0>();
  }
};

template<int Y, int B_V_IB>
using LinearVelocityUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                          ElementVector<>,
                                          ElementVector<Element<Vec3,B_V_IB>>,
                                          MeasLvel>;

template<int B_V_IB>
class LinearVelocityUpdate: public LinearVelocityUpdateBase<0,B_V_IB>{
 public:
  typedef LinearVelocityUpdateBase<0,B_V_IB> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  LinearVelocityUpdate(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = cur.template Get<B_V_IB>() - meas_->GetMeas();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    this->template SetJacCur<0,B_V_IB>(J,cur,Mat3::Identity());
    return 0;
  }
};

} // namespace tsif
