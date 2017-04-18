#ifndef TSIF_POSITION_FINDIF_H_
#define TSIF_POSITION_FINDIF_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int OUT_POS, int STA_POS, int STA_VEL, int STA_ATT>
using PositionFindifBase = Residual<ElementVector<Element<Vec3,OUT_POS>>,
                                 ElementVector<Element<Vec3,STA_POS>,Element<Vec3,STA_VEL>,Element<Quat,STA_ATT>>,
                                 ElementVector<Element<Vec3,STA_POS>>,
                                 MeasEmpty>;

template<int OUT_POS, int STA_POS, int STA_VEL, int STA_ATT>
class PositionFindif: public PositionFindifBase<OUT_POS,STA_POS,STA_VEL,STA_ATT>{
 public:
  typedef PositionFindifBase<OUT_POS,STA_POS,STA_VEL,STA_ATT> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  PositionFindif(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_POS>() = cur.template Get<STA_POS>() - pre.template Get<STA_POS>()
      - dt_*pre.template Get<STA_ATT>().toRotationMatrix()*pre.template Get<STA_VEL>();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_POS)) = -Mat3::Identity();
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_VEL)) = -pre.template Get<STA_ATT>().toRotationMatrix()*dt_;
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_ATT)) = SSM(dt_*pre.template Get<STA_ATT>().toRotationMatrix()*pre.template Get<STA_VEL>());
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_POS),cur.Start(STA_POS)) = Mat3::Identity();
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_POSITION_FINDIF_H_
