#ifndef TSIF_POSE_UPDATE_H_
#define TSIF_POSE_UPDATE_H_

#include "attitude_update.h"
#include "position_update.h"
#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

class MeasPose: public ElementVector<Element<Vec3,0>,Element<Quat,1>>{
 public:
  MeasPose(): ElementVector<Element<Vec3,0>,Element<Quat,1>>(Vec3(0,0,0),Quat(1,0,0,0)){}
  MeasPose(const Vec3& pos, const Quat& att): ElementVector<Element<Vec3,0>,Element<Quat,1>>(pos,att){}
  const Vec3& GetPos() const{
    return Get<0>();
  }
  Vec3& GetPos(){
    return Get<0>();
  }
  const Quat& GetAtt() const{
    return Get<1>();
  }
  Quat& GetAtt(){
    return Get<1>();
  }
};

template<int OUT_POS, int OUT_ATT, int STA_IrIB, int STA_qIB, int STA_IrIJ, int STA_qIJ, int STA_BrBV, int STA_qBV>
using PoseUpdateBase = Residual<ElementVector<Element<Vec3,OUT_POS>,Element<Vec3,OUT_ATT>>,
                                ElementVector<>,
                                ElementVector<Element<Vec3,STA_IrIB>,
                                              Element<Quat,STA_qIB>,
                                              Element<Vec3,STA_IrIJ>,
                                              Element<Quat,STA_qIJ>,
                                              Element<Vec3,STA_BrBV>,
                                              Element<Quat,STA_qBV>>,
                                MeasPose>;

template<int OUT_POS, int OUT_ATT, int STA_IrIB, int STA_qIB, int STA_IrIJ, int STA_qIJ, int STA_BrBV, int STA_qBV>
class PoseUpdate: public PoseUpdateBase<OUT_POS,OUT_ATT,STA_IrIB,STA_qIB,STA_IrIJ,STA_qIJ,STA_BrBV,STA_qBV>{
 public:
  typedef PoseUpdateBase<OUT_POS,OUT_ATT,STA_IrIB,STA_qIB,STA_IrIJ,STA_qIJ,STA_BrBV,STA_qBV> Base;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  PositionUpdate<OUT_POS,STA_IrIB,STA_qIB,STA_IrIJ,STA_qIJ,STA_BrBV> posUpd_;
  AttitudeUpdate<OUT_ATT,STA_qIB,STA_qIJ,STA_qBV> attUpd_;
  PoseUpdate(): Base(false,false,false){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    posUpd_.EvalRes(out,pre,cur);
    attUpd_.EvalRes(out,pre,cur);
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    posUpd_.JacCur(J.block(Output::Start(OUT_POS),0,3,cur.Dim()),pre,cur);
    attUpd_.JacCur(J.block(Output::Start(OUT_ATT),0,3,cur.Dim()),pre,cur);
    return 0;
  }
};

} // namespace tsif

#endif  // TSIF_POSE_UPDATE_H_
