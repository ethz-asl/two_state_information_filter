#ifndef TSIF_POSITION_UPDATE_H_
#define TSIF_POSITION_UPDATE_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasPos: public ElementVector<Element<Vec3,0>>{
 public:
  MeasPos(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasPos(const Vec3& pos): ElementVector<Element<Vec3,0>>(pos){}
  const Vec3& GetPos() const{
    return Get<0>();
  }
  Vec3& GetPos(){
    return Get<0>();
  }
};

template<int OUT_POS, int STA_IrIB, int STA_qIB, int STA_IrIJ, int STA_qIJ, int STA_BrBV>
using PositionUpdateBase = Residual<ElementVector<Element<Vec3,OUT_POS>>,
                                    ElementVector<>,
                                    ElementVector<Element<Vec3,STA_IrIB>,
                                                  Element<Quat,STA_qIB>,
                                                  Element<Vec3,STA_IrIJ>,
                                                  Element<Quat,STA_qIJ>,
                                                  Element<Vec3,STA_BrBV>>,
                                    MeasPos>;

template<int OUT_POS, int STA_IrIB, int STA_qIB, int STA_IrIJ, int STA_qIJ, int STA_BrBV>
class PositionUpdate: public PositionUpdateBase<OUT_POS,STA_IrIB,STA_qIB,STA_IrIJ,STA_qIJ,STA_BrBV>{
 public:
  typedef PositionUpdateBase<OUT_POS,STA_IrIB,STA_qIB,STA_IrIJ,STA_qIJ,STA_BrBV> Base;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  PositionUpdate(): Base(false,false,false){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_POS>() = cur.template Get<STA_qIJ>().inverse().toRotationMatrix()*(
        cur.template Get<STA_IrIB>() - cur.template Get<STA_IrIJ>()
       + cur.template Get<STA_qIB>().toRotationMatrix()*cur.template Get<STA_BrBV>()) - meas_->GetPos();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 pos = cur.template Get<STA_qIJ>().inverse().toRotationMatrix()*(
        cur.template Get<STA_IrIB>() - cur.template Get<STA_IrIJ>()
       + cur.template Get<STA_qIB>().toRotationMatrix()*cur.template Get<STA_BrBV>());
    const Mat3 mJI =  cur.template Get<STA_qIJ>().inverse().toRotationMatrix();
    const Mat3 mIB =  cur.template Get<STA_qIB>().toRotationMatrix();
    this->template SetJacCur<OUT_POS,STA_IrIB>(J,cur,mJI);
    this->template SetJacCur<OUT_POS,STA_qIB>(J,cur,-mJI*SSM(mIB*cur.template Get<STA_BrBV>()));
    this->template SetJacCur<OUT_POS,STA_IrIJ>(J,cur,-mJI);
    this->template SetJacCur<OUT_POS,STA_qIJ>(J,cur,SSM(pos)*mJI);
    this->template SetJacCur<OUT_POS,STA_BrBV>(J,cur,mJI*mIB);
    return 0;
  }
};

} // namespace tsif

#endif  // TSIF_POSITION_UPDATE_H_
