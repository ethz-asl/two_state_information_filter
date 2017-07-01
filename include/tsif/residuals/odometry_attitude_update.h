#pragma once

#include "tsif/residual.h"

namespace tsif{

//measurement class for for odom to base orientation
class MeasAttOdom: public ElementVector<Element<Quat,0>>{
 public:
  //constructors with and without arguments
  MeasAttOdom(): ElementVector<Element<Quat,0>>(Quat(1,0,0,0)){}
  MeasAttOdom(const Quat& q_IB): ElementVector<Element<Quat,0>>(q_IB){}
  //functions to get  the measurement
  const Quat& GetAtt() const{
    return Get<0>();
  }
  Quat& GetAtt(){
    return Get<0>();
  }
};

//base class type for convenience
template<int Y,//innovation
         int PHI_JB,//map to base orientation
         int PHI_JI>//map to odom orientation
using OdometryAttitudeUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                              ElementVector<>,
                                              ElementVector<Element<Quat,PHI_JB>,
                                                            Element<Quat,PHI_JI>>,
                                              MeasAttMapCentric>;

template<int PHI_JB,int PHI_JI>
class OdometryAttitudeUpdate: public OdometryAttitudeUpdateBase<0,PHI_JB,PHI_JI>{
 public:
  //types from base class
  typedef OdometryAttitudeUpdateBase<0,PHI_JB,PHI_JI> Base;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  //measurement pointer from base class
  using Base::meas_;
  //constructor
  OdometryAttitudeUpdate(): Base(false,false,false){}
  //function to evaluate the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = Log(cur.template Get<PHI_JI>()*meas_->GetAtt()*cur.template Get<PHI_JB>().inverse());
    return 0;
  }
  //function evaluating the innovation jacobian wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    //has dimension 0
    return 0;
  }
  //function evaluating the innovation jacobian wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 y = Log(cur.template Get<PHI_JI>()*meas_->GetAtt()*cur.template Get<PHI_JB>().inverse());
    const Mat3 C_JB =  cur.template Get<PHI_JI>().toRotationMatrix()*meas_->GetAtt().toRotationMatrix();
    const Mat3 G_inv = GammaMat(y).inverse();
    this->template SetJacCur<0,PHI_JB>(J,cur,-G_inv*C_JB*cur.template Get<PHI_JB>().inverse().toRotationMatrix());
    this->template SetJacCur<0,PHI_JI>(J,cur,G_inv);
    return 0;
  }
};

} // namespace tsif
