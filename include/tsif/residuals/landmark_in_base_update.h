#pragma once

//base class
#include "tsif/residual.h"

namespace tsif {

//measurement class to hold base landmark vector in base frame
class LandmarkInBaseMeasurement : public ElementVector<Element<Vec3,0>>{
 public:
  //constructors with and without arguments
  LandmarkInBaseMeasurement(): ElementVector<Element<Vec3,0>>(Vec3(.0,.0,.0)) {}
  LandmarkInBaseMeasurement(const Vec3& B_Ns):ElementVector<Element<Vec3,0>>(B_Ns){}
  //functions to measurement
  const Vec3& GetMeasurement() const {return Get<0>();}
  Vec3& GetMeasurement() {return Get<0>();}
};

//shortcut to make the templating more convenient
template<int Y,//innovation
         int B_P>//landmark position in base frame
using LandmarkInBaseUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,//innovation
                                          ElementVector<>,//previous state is empty
                                          ElementVector<Element<Vec3,B_P>>,//current state
                                          LandmarkInBaseMeasurement>;//measurement

//update residual comparing a measured landmark position to a landmark in the filter state
template<int B_P>
class LandmarkInBaseUpdate : public LandmarkInBaseUpdateBase<0,B_P> {
 public:
  //typedef and types from base for convenience
  //base type
  typedef LandmarkInBaseUpdateBase<0,B_P> Base;
  //innovation vector type
  using typename Base::Output;
  //previous state vector type
  using typename Base::Previous;
  //current state vector type
  using typename Base::Current;
  //measurement pointer from base
  using Base::meas_;
  //constructor setting characteristics of the residual and default weights
  LandmarkInBaseUpdate(): Base (false,//splittable
                                false,//mergeable
                                false){};//mandatory

  //function evaluating the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
   //compare measured contact points to filter state
		 out. template Get<0>() = meas_->GetMeasurement()-cur.template Get<B_P>();
   return 0;
  }
  //function evaluating the jacobian of the the innovation wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //has dimension 0
   return 0;    
  }
  //function evaluating the jacobian of the the innovation wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //set jacobian zero
   J.setZero();
   //rotation matrix from odom to base frame from orientation state   
   this->template SetJacCur<0,B_P>(J,cur,-Mat3::Identity());
   return 0;    
  }
};

} // namespace tsif

