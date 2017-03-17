#ifndef TSIF_GYROSCOPE_UPDATE_H_
#define TSIF_GYROSCOPE_UPDATE_H_

#include "tsif/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasGyr: public ElementVector<Element<Vec3,0>>{
 public:
  MeasGyr(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasGyr(const Vec3& acc): ElementVector<Element<Vec3,0>>(acc){}
  const Vec3& GetGyr() const{
    return Get<0>();
  }
  Vec3& GetGyr(){
    return Get<0>();
  }
};

template<int OUT_ROR, int STA_ROR, int STA_GYB>
using GyroscopeUpdateBase = Residual<ElementVector<Element<Vec3,OUT_ROR>>,
                                     ElementVector<>,
                                     ElementVector<Element<Vec3,STA_ROR>,Element<Vec3,STA_GYB>>,
                                     MeasGyr>;

template<int OUT_ROR, int STA_ROR, int STA_GYB>
class GyroscopeUpdate: public GyroscopeUpdateBase<OUT_ROR,STA_ROR,STA_GYB>{
 public:
  typedef GyroscopeUpdateBase<OUT_ROR,STA_ROR,STA_GYB> Base;
  using Base::dt_;
  using Base::w_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  GyroscopeUpdate(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_ROR>() = cur.template Get<STA_ROR>() + cur.template Get<STA_GYB>() - meas_->GetGyr();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_ROR),cur.Start(STA_ROR)) = Mat3::Identity();
    J.block<3,3>(Output::Start(OUT_ROR),cur.Start(STA_GYB)) = Mat3::Identity();
    return 0;
  }
  double GetWeight(){
    return w_*sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_GYROSCOPE_UPDATE_H_
