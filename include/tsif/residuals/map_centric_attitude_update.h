#pragma once

#include "tsif/residual.h"
#include "tsif/utils/common.h"

namespace tsif{

//measurement class for for measurement frame to map orientation
class MeasAttMapCentric: public ElementVector<Element<Quat,0>>{
 public:
  //constructors with and without arguments
  MeasAttMapCentric(): ElementVector<Element<Quat,0>>(Quat(1,0,0,0)){}
  MeasAttMapCentric(const Quat& q_JV): ElementVector<Element<Quat,0>>(q_JV){}
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
         int PHI_BV>//mbase to measurement frame orientation
using MapCentricAttitudeUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                              ElementVector<>,
                                              ElementVector<Element<Quat,PHI_JB>,
                                                            Element<Quat,PHI_BV>>,
                                              MeasAttMapCentric>;

template<int PHI_JB,int PHI_BV>
class MapCentricAttitudeUpdate: public MapCentricAttitudeUpdateBase<0,PHI_JB,PHI_BV>{
 public:
  //types from base class
  typedef MapCentricAttitudeUpdateBase<0,PHI_JB,PHI_BV> Base;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  //measurement pointer from base class
  using Base::meas_;
  using Base::w_;
  double huber_threshold_;
  //constructor
  MapCentricAttitudeUpdate(): Base(false,false,false), huber_threshold_(1000.){}
  //function to evaluate the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = Log(cur.template Get<PHI_JB>()*cur.template Get<PHI_BV>()*meas_->GetAtt().inverse());
    return 0;
  }
  //function evaluating the innovation jacobian wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    //has dimension 0
    return 0;
  }
  //function evaluating the innovation jacobian wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 y = Log(cur.template Get<PHI_JB>()*cur.template Get<PHI_BV>()*meas_->GetAtt().inverse());
    const Mat3 C_JB =  cur.template Get<PHI_JB>().toRotationMatrix();
    const Mat3 G_inv = GammaMat(y).inverse();
    this->template SetJacCur<0,PHI_JB>(J,cur,G_inv);
    this->template SetJacCur<0,PHI_BV>(J,cur,G_inv*C_JB);
    return 0;
  }
  virtual void AddNoise(typename Output::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    //compute weight using huber loss function
    double weight = w_;
    const double norm = out.template Get<0>().norm();
    if(norm > huber_threshold_) weight *= sqrt(2 * huber_threshold_ * (norm - 0.5 * huber_threshold_)/(norm*norm));
    //scale the innovation and jacobians
    out.Scale(weight);
    J_pre *= weight;
    J_cur *= weight;
  }
};

} // namespace tsif
