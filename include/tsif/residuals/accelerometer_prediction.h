#ifndef TSIF_ACCELEROMETER_PREDICTION_H_
#define TSIF_ACCELEROMETER_PREDICTION_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasAcc: public ElementVector<Element<Vec3,0>>{
 public:
  MeasAcc(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasAcc(const Vec3& acc): ElementVector<Element<Vec3,0>>(acc){}
  const Vec3& GetAcc() const{
    return Get<0>();
  }
  Vec3& GetAcc(){
    return Get<0>();
  }
};

template<int OUT_VEL, int STA_VEL, int STA_ATT, int STA_ROR, int STA_ACB>
using AccelerometerPredictionBase = Residual<ElementVector<Element<Vec3,OUT_VEL>>,
                                           ElementVector<Element<Vec3,STA_VEL>,Element<Quat,STA_ATT>,Element<Vec3,STA_ROR>,Element<Vec3,STA_ACB>>,
                                           ElementVector<Element<Vec3,STA_VEL>>,
                                           MeasAcc>;

template<int OUT_VEL, int STA_VEL, int STA_ATT, int STA_ROR, int STA_ACB>
class AccelerometerPrediction: public AccelerometerPredictionBase<OUT_VEL,STA_VEL,STA_ATT,STA_ROR,STA_ACB>{
 public:
  typedef AccelerometerPredictionBase<OUT_VEL,STA_VEL,STA_ATT,STA_ROR,STA_ACB> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  const Vec3 g_;
  AccelerometerPrediction(bool isSplitable = true,bool isMergeable = true,bool isMandatory = true): Base(isSplitable,isMergeable,isMandatory), g_(0,0,-9.81){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_VEL>() = cur.template Get<STA_VEL>()
        - (Mat3::Identity() - SSM(dt_*pre.template Get<STA_ROR>()))*pre.template Get<STA_VEL>()
        - dt_*(meas_->GetAcc()-pre.template Get<STA_ACB>()+pre.template Get<STA_ATT>().inverse().toRotationMatrix()*g_);
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    this->template SetJacPre<OUT_VEL, STA_VEL>(J, pre, -(Mat3::Identity() - SSM(dt_*pre.template Get<STA_ROR>())));
    this->template SetJacPre<OUT_VEL, STA_ATT>(J, pre, -pre.template Get<STA_ATT>().inverse().toRotationMatrix()*SSM(g_*dt_));
    this->template SetJacPre<OUT_VEL, STA_ROR>(J, pre, -SSM(dt_*pre.template Get<STA_VEL>()));
    this->template SetJacPre<OUT_VEL, STA_ACB>(J, pre, dt_*Mat3::Identity());
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    this->template SetJacCur<OUT_VEL, STA_VEL>(J, cur, Mat3::Identity());
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_ACCELEROMETER_PREDICTION_H_
