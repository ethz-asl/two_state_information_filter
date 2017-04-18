#ifndef TSIF_ATTITUDE_UPDATE_H_
#define TSIF_ATTITUDE_UPDATE_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasAtt: public ElementVector<Element<Quat,0>>{
 public:
  MeasAtt(): ElementVector<Element<Quat,0>>(Quat(1,0,0,0)){}
  MeasAtt(const Quat& att): ElementVector<Element<Quat,0>>(att){}
  const Quat& GetAtt() const{
    return Get<0>();
  }
  Quat& GetAtt(){
    return Get<0>();
  }
};

template<int OUT_ATT, int STA_qIB, int STA_qIJ, int STA_qBV>
using AttitudeUpdateBase = Residual<ElementVector<Element<Vec3,OUT_ATT>>,
                                    ElementVector<>,
                                    ElementVector<Element<Quat,STA_qIB>,
                                                  Element<Quat,STA_qIJ>,
                                                  Element<Quat,STA_qBV>>,
                                    MeasAtt>;

template<int OUT_ATT, int STA_qIB, int STA_qIJ, int STA_qBV>
class AttitudeUpdate: public AttitudeUpdateBase<OUT_ATT,STA_qIB,STA_qIJ,STA_qBV>{
 public:
  typedef AttitudeUpdateBase<OUT_ATT,STA_qIB,STA_qIJ,STA_qBV> Base;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  AttitudeUpdate(): Base(false,false,false){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_ATT>() = Log(cur.template Get<STA_qIJ>().inverse()*
                                      cur.template Get<STA_qIB>()*cur.template Get<STA_qBV>()*
                                      meas_->GetAtt().inverse());
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 attErr = Log(cur.template Get<STA_qIJ>().inverse()*
                                      cur.template Get<STA_qIB>()*cur.template Get<STA_qBV>()*
                                      meas_->GetAtt().inverse());
    const Mat3 mJI =  cur.template Get<STA_qIJ>().inverse().toRotationMatrix();
    const Mat3 mIB =  cur.template Get<STA_qIB>().toRotationMatrix();
    const Mat3 GI = GammaMat(attErr).inverse();
    this->template SetJacCur<OUT_ATT,STA_qIB>(J,cur,GI*mJI);
    this->template SetJacCur<OUT_ATT,STA_qIJ>(J,cur,-GI*mJI);
    this->template SetJacCur<OUT_ATT,STA_qBV>(J,cur,GI*mJI*mIB);
    return 0;
  }
};

} // namespace tsif

#endif  // TSIF_ATTITUDE_UPDATE_H_
