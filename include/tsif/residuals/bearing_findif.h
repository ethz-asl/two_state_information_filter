#ifndef TSIF_BEARING_FINDIF_H_
#define TSIF_BEARING_FINDIF_H_

#include "tsif/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int OUT_BEA, int STA_BEA, int STA_DIS, int STA_VEL, int STA_ROR, int N>
using BearingFindifBase = Residual<ElementVector<Element<std::array<Vec<2>,N>,OUT_BEA>>,
                                   ElementVector<Element<std::array<UnitVector,N>,STA_BEA>,Element<std::array<double,N>,STA_DIS>,Element<Vec3,STA_VEL>,Element<Vec3,STA_ROR>>,
                                   ElementVector<Element<std::array<UnitVector,N>,STA_BEA>>,
                                   MeasEmpty>;

template<int OUT_BEA, int STA_BEA, int STA_DIS, int STA_VEL, int STA_ROR, int N>
class BearingFindif: public BearingFindifBase<OUT_BEA,STA_BEA,STA_DIS,STA_VEL,STA_ROR,N>{
 public:
  typedef BearingFindifBase<OUT_BEA,STA_BEA,STA_DIS,STA_VEL,STA_ROR,N> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  BearingFindif(): Base(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    UnitVector n_predicted;
    const Vec3& ror = pre.template Get<STA_ROR>();
    const Vec3& vel = pre.template Get<STA_VEL>();
    for(int i=0;i<N;i++){
      const UnitVector& bea = pre.template Get<STA_BEA>()[i];
      const Vec3 beaVec = bea.GetVec();
      const Mat<3,2> beaN = bea.GetN();
      const double& invDis = pre.template Get<STA_DIS>()[i];
      const Vec<2> dn = beaN.transpose()*(ror + invDis * beaVec.cross(vel));
      bea.Boxplus(dn*dt_,n_predicted);
      n_predicted.Boxminus(cur.template Get<STA_BEA>()[i],out.template Get<OUT_BEA>()[i]); // CAREFUL: DO NOT INVERSE BOXMINUS ORDER (CHAIN-RULE NOT VALID)
    }
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    UnitVector n_predicted;
    const Vec3& ror = pre.template Get<STA_ROR>();
    const Vec3& vel = pre.template Get<STA_VEL>();
    for(int i=0;i<N;i++){
      const UnitVector& bea = pre.template Get<STA_BEA>()[i];
      const Vec3 beaVec = bea.GetVec();
      const Mat<3,2> beaN = bea.GetN();
      const Mat<3,2> beaM = bea.GetM();
      const double& invDis = pre.template Get<STA_DIS>()[i];
      const Vec<2> dn = beaN.transpose()*(ror + invDis * beaVec.cross(vel));
      bea.Boxplus(dn*dt_,n_predicted);
      Mat<2> Jsub_1,Jsub_2a, Jsub_2b;
      bea.BoxplusJacInp(dn*dt_,Jsub_2a);
      bea.BoxplusJacVec(dn*dt_,Jsub_2b);
      n_predicted.BoxminusJacInp(cur.template Get<STA_BEA>()[i],Jsub_1);
      J.block<2,3>(Output::Start(OUT_BEA)+2*i,pre.Start(STA_ROR)) = dt_*Jsub_1*Jsub_2b*beaN.transpose();
      J.block<2,3>(Output::Start(OUT_BEA)+2*i,pre.Start(STA_VEL)) = dt_*invDis*Jsub_1*Jsub_2b*beaN.transpose()*SSM(beaVec);
      J.block<2,1>(Output::Start(OUT_BEA)+2*i,pre.Start(STA_DIS)+i) = dt_*Jsub_1*Jsub_2b*beaN.transpose()*beaVec.cross(vel);
      J.block<2,2>(Output::Start(OUT_BEA)+2*i,pre.Start(STA_BEA)+2*i) = dt_*Jsub_1*Jsub_2b*(-invDis*beaN.transpose()*SSM(vel)*beaM
          +beaN.transpose()*SSM(ror + invDis * beaVec.cross(vel))*beaN) + Jsub_1*Jsub_2a;
    }
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    UnitVector n_predicted;
    const Vec3& ror = pre.template Get<STA_ROR>();
    const Vec3& vel = pre.template Get<STA_VEL>();
    for(int i=0;i<N;i++){
      const UnitVector& bea = pre.template Get<STA_BEA>()[i];
      const Vec3 beaVec = bea.GetVec();
      const Mat<3,2> beaN = bea.GetN();
      const double& invDis = pre.template Get<STA_DIS>()[i];
      const Vec<2> dn = beaN.transpose()*(ror + invDis * beaVec.cross(vel));
      bea.Boxplus(dn*dt_,n_predicted);
      Mat<2> Jsub_1;
      n_predicted.BoxminusJacRef(cur.template Get<STA_BEA>()[i],Jsub_1);
      J.block<2,2>(Output::Start(OUT_BEA)+2*i,cur.Start(STA_BEA)+2*i) = Jsub_1;
    }
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_BEARING_FINDIF_H_
