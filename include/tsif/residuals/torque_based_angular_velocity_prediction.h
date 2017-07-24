#pragma once

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{
class MeasExternalForces: public ElementVector<Element<std::array<Vec3,4>,0>,
                                               Element<Vec<12>,1>>{
 public:
  MeasExternalForces(){
    for (int i=0; i<Get<0>().size();i++) Get<0>()[i]=Vec3::Zero();
    Get<1>()=Vec<12>::Zero();
  }
  MeasExternalForces(const Vec3& B_f_rf, const Vec3& B_f_lf, const Vec3& B_f_lh, const Vec3& B_f_rh, const Vec<12>& theta){
    Get<0>()[0]=B_f_rf;
    Get<0>()[1]=B_f_lf;
    Get<0>()[2]=B_f_lh;
    Get<0>()[3]=B_f_rh;
    Get<1>()=theta;
  }
  const Vec3& GetForce(int leg_index) const{
    return Get<0>()[leg_index];
  }
  Vec3& GetForce(int leg_index) {
    return Get<0>()[leg_index];
  }
  const Vec<12>& GetTorques() const{
    return Get<1>();
  }
  Vec<12>& GetTorques(){
    return Get<1>();
  }
};

template<int Y, int B_OMEGA_IB, int B_COM_B>
using TorqueBasedAngularVelocityPredictionBase = Residual<ElementVector<Element<Vec3,Y>>,
                                                          ElementVector<Element<Vec3,B_COM_B>>,
                                                          ElementVector<Element<Vec3,B_OMEGA_IB>>,
                                                          MeasExternalForces>;

template<int B_OMEGA_IB, int B_COM_B>
class TorqueBasedAngularVelocityPrediction: public TorqueBasedAngularVelocityPredictionBase<0, B_OMEGA_IB, B_COM_B>{
 public:
  typedef TorqueBasedAngularVelocityPredictionBase<0, B_OMEGA_IB, B_COM_B> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  TorqueBasedAngularVelocityPrediction(): Base(true,true,false){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif
