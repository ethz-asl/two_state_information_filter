#pragma once

#include "tsif/residual.h"

namespace tsif{

//measurement class for the odometry position estimate
class MeasPosOdom: public ElementVector<Element<Vec3,0>>{
 public:
  //constructors with and without arguments
  MeasPosOdom(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasPosOdom(const Vec3& I_r_IB): ElementVector<Element<Vec3,0>>(I_r_IB){}
  //functions to get the position measurement
  const Vec3& GetPos() const{
    return Get<0>();
  }
  Vec3& GetPos(){
    return Get<0>();
  }
};

//base class type for convenience
template<int Y,//innovation
         int J_R_JB,//map to base position in map frame
         int J_R_JI,//map to odom position in map frame
         int PHI_JI>//map to odom orientation
using OdometryPositionUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                              ElementVector<>,
                                              ElementVector<Element<Vec3,J_R_JB>,
                                                            Element<Vec3,J_R_JI>,
                                                            Element<Quat,PHI_JI>>,
                                              MeasPosOdom>;

//residual class implementing an update of the odom frame position in map frame
template<int J_R_JB, int J_R_JI, int PHI_JI>
class OdometryPositionUpdate: public OdometryPositionUpdateBase<0,J_R_JB,J_R_JI,PHI_JI>{
 public:
  //types from base class 
  typedef OdometryPositionUpdateBase<0,J_R_JB,J_R_JI,PHI_JI> Base;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  //measurement pointer from base class
  using Base::meas_;
  //constructor
  OdometryPositionUpdate(): Base(false,false,false){}
  //function evaluating the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = cur.template Get<J_R_JI>() - cur.template Get<J_R_JB>() + cur.template Get<PHI_JI>().toRotationMatrix()*meas_->GetPos();
    return 0;
  }
  //function evaluating the innovation jacobian wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    //has dimension zero
    return 0;
  }
  //function evaluating the innovation jacobian wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    this->template SetJacCur<0,J_R_JI>(J,cur,Mat3::Identity());
    this->template SetJacCur<0,J_R_JB>(J,cur,-Mat3::Identity());
    this->template SetJacCur<0,PHI_JI>(J,cur,-SSM(cur.template Get<PHI_JI>()*meas_->GetPos()));
    return 0;
  }
};

} // namespace tsif
