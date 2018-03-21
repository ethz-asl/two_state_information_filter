#ifndef TSIF_DISTANCE_FINDIF_H_
#define TSIF_DISTANCE_FINDIF_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int OUT_DIS, int STA_BEA, int STA_DIS, int STA_VEL, int STA_ROR, int STA_VEP, int STA_VEA, int N>
using DistanceFindifBase = Residual<ElementVector<Element<std::array<Vec<1>,N>,OUT_DIS>>,
                                    ElementVector<Element<std::array<UnitVector,N>,STA_BEA>,
                                                  Element<std::array<double,N>,STA_DIS>,
                                                  Element<Vec3,STA_VEL>,
                                                  Element<Vec3,STA_ROR>,
                                                  Element<Vec3,STA_VEP>,
                                                  Element<Quat,STA_VEA>>,
                                    ElementVector<Element<std::array<double,N>,STA_DIS>>,
                                    MeasEmpty>;

template<int OUT_DIS, int STA_BEA, int STA_DIS, int STA_VEL, int STA_ROR, int STA_VEP, int STA_VEA, int N>
class DistanceFindif: public DistanceFindifBase<OUT_DIS,STA_BEA,STA_DIS,STA_VEL,STA_ROR,STA_VEP,STA_VEA,N>{
 public:
  typedef DistanceFindifBase<OUT_DIS,STA_BEA,STA_DIS,STA_VEL,STA_ROR,STA_VEP,STA_VEA,N> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  DistanceFindif(): Base(true,true,true), vep_not_fixed_(STA_VEP>=0), vea_not_fixed_(STA_VEA>=0){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    const Mat3 C_VI = pre.template Get<STA_VEA>().toRotationMatrix();
    const Vec3 vel = C_VI*(pre.template Get<STA_VEL>() + pre.template Get<STA_ROR>().cross(pre.template Get<STA_VEP>()));
    for(int i=0;i<N;i++){
      const UnitVector& bea = pre.template Get<STA_BEA>()[i];
      const Vec3 beaVec = bea.GetVec();
      const double& invDis = pre.template Get<STA_DIS>()[i];
      out.template Get<OUT_DIS>()[i](0) = invDis + dt_*beaVec.dot(vel)*invDis*invDis - cur.template Get<STA_DIS>()[i];
    }
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    const Mat3 C_VI = pre.template Get<STA_VEA>().toRotationMatrix();
    const Vec3 vel = C_VI*(pre.template Get<STA_VEL>() + pre.template Get<STA_ROR>().cross(pre.template Get<STA_VEP>()));
    for(int i=0;i<N;i++){
      const UnitVector& bea = pre.template Get<STA_BEA>()[i];
      const Vec3 beaVec = bea.GetVec();
      const double& invDis = pre.template Get<STA_DIS>()[i];
      const Mat<1,3> J_vel = dt_*beaVec.transpose()*invDis*invDis;
      J.block<1,3>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_VEL)) = J_vel*C_VI;
      J.block<1,3>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_ROR)) = -J_vel*C_VI*SSM(pre.template Get<STA_VEP>());
      J.block<1,1>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_DIS)+1*i) = Mat<1>::Identity()+dt_*2*beaVec.transpose()*vel*invDis;
      J.block<1,2>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_BEA)+2*i) = dt_*vel.transpose()*invDis*invDis*bea.GetM();
      if (vep_not_fixed_) J.block<1,3>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_VEP)) = J_vel*C_VI*SSM(pre.template Get<STA_ROR>());
      if (vea_not_fixed_) J.block<1,3>(Output::Start(OUT_DIS)+1*i,pre.Start(STA_VEA)) = - J_vel*SSM(vel);
    }
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.setZero();
    for(int i=0;i<N;i++){
      J.block<1,1>(Output::Start(OUT_DIS)+1*i,cur.Start(STA_DIS)+1*i) = -Mat<1>::Identity();
    }
    return 0;
  }
  double GetWeight(){
    return w_/sqrt(dt_);
  }
 protected:
  const bool vep_not_fixed_;
  const bool vea_not_fixed_;
};

} // namespace tsif

#endif  // TSIF_DISTANCE_FINDIF_H_
