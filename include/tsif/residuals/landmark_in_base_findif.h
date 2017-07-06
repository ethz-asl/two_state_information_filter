#pragma once

//base class
#include "tsif/residual.h"

namespace tsif {

//measurement class to hold a bool indicating that the landmark is active (so that the weighting can be adjusted accordingly)
class LandmarkInBaseFlagMeasurement : public ElementVector<Element<bool,0>>{
 public:
  LandmarkInBaseFlagMeasurement() : ElementVector<Element<bool,0>>(false){};
  LandmarkInBaseFlagMeasurement(const bool& b) : ElementVector<Element<bool,0>>(b){};
  const bool& GetMeasurement() const {return Get<0>();};
  bool& GetMeasurement() {return Get<0>();};
};

//shortcut to make the templating more convenient
template<int Y,//contact point innovation
         int B_V_IB,//linear velocity base to odom in base frame
         int B_OMEGA_IB,//angular velocity base to odom in base
         int B_P>//landmark in base frame
using LandmarkInBaseFindifBase = Residual<ElementVector<Element<Vec3,Y>>,//innovation
                                    ElementVector<Element<Vec3,B_V_IB>,Element<Vec3,B_OMEGA_IB>,Element<Vec3,B_P>>,//previous state
                                    ElementVector<Element<Vec3,B_P>>,//current state
                                    LandmarkInBaseFlagMeasurement>;//measurement

//prediction residual comparing current landmark position with 
template<int B_V_IB,int B_OMEGA_IB,int B_P>
class LandmarkInBaseFindif : public LandmarkInBaseFindifBase<0,B_V_IB,B_OMEGA_IB,B_P> {
 public:
  //typedef and types from base for convenience
  //base type
  typedef LandmarkInBaseFindifBase<0,B_V_IB,B_OMEGA_IB,B_P> Base;
  //innovation vector type
  using typename Base::Output;
  //previous state vector type
  using typename Base::Previous;
  //current state vector type
  using typename Base::Current;
  //time increment from base
  using Base::dt_;
  //measurement pointer from base
  using Base::meas_;
  //weights
  double w_active_;
  double w_passive_;
  //constructor setting characteristics of the residual and default weights
  LandmarkInBaseFindif(): Base (true,//splittable
                                true,//mergeable
                                true),//mandatory
                                w_active_(1.), w_passive_(0.000001) {};
  //function evaluating the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
   //compare measured contact points to filter state
		 out. template Get<0>()=cur.template Get<B_P>()
                          -(Mat3::Identity()+SSM(dt_*pre.template Get<B_OMEGA_IB>()))*pre.template Get<B_P>()
                          +dt_*pre.template Get<B_V_IB>();
   return 0;
  }
  //function evaluating the jacobian of the the innovation wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //set jacobian zero
   J.setZero();
   //set nonzero blocks of the jacobian
   //derivative of innovation wrt linear velocity
   this->template SetJacPre<0,B_V_IB>(J,pre,dt_*Mat3::Identity());
   //derivative of innovation wrt angular velocity
   this->template SetJacPre<0,B_OMEGA_IB>(J,pre,SSM(dt_*pre.template Get<B_P>()));
   //derivative of innovation wrt landmark
   this->template SetJacPre<0,B_P>(J,pre,-SSM(dt_*pre.template Get<B_OMEGA_IB>())-Mat3::Identity());
   return 0;    
  }
  //function evaluating the jacobian of the the innovation wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //set jacobian zero
   J.setZero();
   //set nonzero blocks of the jacobian
   //derivative of innovation wrt linear velocity
   this->template SetJacCur<0,B_P>(J,cur,Mat3::Identity());
   return 0;    
  }
  //function returning the weighting factor used in the TSIF optimization
  double GetWeight(){
    if (meas_->GetMeasurement()) return w_active_/sqrt(dt_);
    else return w_passive_/sqrt(dt_);
  }
};

} // namespace tsif

