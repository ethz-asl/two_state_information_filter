#pragma once

#include "tsif/residual.h"
#include "tsif/utils/common.h"

namespace tsif{

//measurement class for the measurement frame position in map frame
class MeasPosMapCentric: public ElementVector<Element<Vec3,0>>{
 public:
  //constructors with and without arguments
  MeasPosMapCentric(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasPosMapCentric(const Vec3& J_r_JV): ElementVector<Element<Vec3,0>>(J_r_JV){}
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
         int PHI_JB,//map to base orientation
         int B_R_BV>//position of measurement frame in base
using MapCentricPositionUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,
                                              ElementVector<>,
                                              ElementVector<Element<Vec3,J_R_JB>,
                                                            Element<Quat,PHI_JB>,
                                                            Element<Vec3,B_R_BV>>,
                                              MeasPosMapCentric>;

//residual class implementing an update of the position in the map frame from a direct measurement
template<int J_R_JB, int PHI_JB, int B_R_BV>
class MapCentricPositionUpdate: public MapCentricPositionUpdateBase<0,J_R_JB,PHI_JB,B_R_BV>{
 public:
  //types from base class 
  typedef MapCentricPositionUpdateBase<0,J_R_JB,PHI_JB,B_R_BV> Base;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  //measurement pointer from base class
  using Base::meas_;
  using Base::w_;
  double huber_threshold_;
  //constructor
  MapCentricPositionUpdate(): Base(false,false,false), huber_threshold_(1000.){}
  //function evaluating the innovation
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<0>() = cur.template Get<J_R_JB>() + cur.template Get<PHI_JB>().toRotationMatrix()*cur.template Get<B_R_BV>()
                            -meas_->GetPos();
    return 0;
  }
  //function evaluating the innovation jacobian wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    //has dimension zero
    return 0;
  }
  //function evaluating the innovation jacobian wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Mat3 C_JB = cur.template Get<PHI_JB>().toRotationMatrix();
    this->template SetJacCur<0,J_R_JB>(J,cur,Mat3::Identity());
    this->template SetJacCur<0,PHI_JB>(J,cur,-SSM(C_JB*cur.template Get<B_R_BV>()));
    this->template SetJacCur<0,B_R_BV>(J,cur,C_JB);
    return 0;
  }
  //function to scale the innovation and jacobians with the appropriate weighting factor
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
