#ifndef TSIF_ATTITUDE_FINDIF_H_
#define TSIF_ATTITUDE_FINDIF_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int OUT_ATT, int STA_ATT, int STA_ROR>
using AttitudeFindifBase = Residual<ElementVector<Element<Vec3,OUT_ATT>>,
                                 ElementVector<Element<Quat,STA_ATT>,Element<Vec3,STA_ROR>>,
                                 ElementVector<Element<Quat,STA_ATT>>,
                                 MeasEmpty>;

template<int OUT_ATT, int STA_ATT, int STA_ROR>
class AttitudeFindif: public AttitudeFindifBase<OUT_ATT,STA_ATT,STA_ROR>{
 public:
  typedef AttitudeFindifBase<OUT_ATT,STA_ATT,STA_ROR> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  AttitudeFindif(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_ATT>() = Boxminus(Boxplus(pre.template Get<STA_ATT>(),dt_*pre.template Get<STA_ROR>()),
                                           cur.template Get<STA_ATT>());
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 err = Boxminus(Boxplus(pre.template Get<STA_ATT>(),dt_*pre.template Get<STA_ROR>()),
                                           cur.template Get<STA_ATT>());
    J.block<3,3>(Output::Start(OUT_ATT),pre.Start(STA_ATT)) =
        GammaMat(err).inverse()*(Exp(dt_*pre.template Get<STA_ROR>())).toRotationMatrix();
    J.block<3,3>(Output::Start(OUT_ATT),pre.Start(STA_ROR)) =
        dt_*GammaMat(err).inverse()*GammaMat(dt_*pre.template Get<STA_ROR>());
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Vec3 err = Boxminus(Boxplus(pre.template Get<STA_ATT>(),dt_*pre.template Get<STA_ROR>()),
                                           cur.template Get<STA_ATT>());
    J.block<3,3>(Output::Start(OUT_ATT),cur.Start(STA_ATT)) = -GammaMat(err).inverse()
                                                              *Boxplus(pre.template Get<STA_ATT>(),dt_*pre.template Get<STA_ROR>()).toRotationMatrix()
                                                              *cur.template Get<STA_ATT>().toRotationMatrix().transpose();
    return 0;
  }
  virtual double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_ATTITUDE_FINDIF_H_
