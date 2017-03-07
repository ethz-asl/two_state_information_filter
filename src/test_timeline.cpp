#include "tsif/filter.h"

using namespace tsif;

class MeasA: public ElementVector<Element<Vec3,0>>{
 public:
  MeasA(): ElementVector<Element<Vec3,0>>(Vec3(0,0,0)){}
  MeasA(const Vec3& acc): ElementVector<Element<Vec3,0>>(acc){}
  const Vec3& GetAcc() const{
    return Get<0>();
  }
  Vec3& GetAcc(){
    return Get<0>();
  }
};

template<int OUT_POS, int STA_POS, int STA_VEL>
using ResidualPosBase = Residual<ElementVector<Element<Vec3,OUT_POS>>,
                                 ElementVector<Element<Vec3,STA_POS>,Element<Vec3,STA_VEL>>,
                                 ElementVector<Element<Vec3,STA_POS>>,
                                 MeasEmpty>;

template<int OUT_POS, int STA_POS, int STA_VEL>
class ResidualPos: public ResidualPosBase<OUT_POS,STA_POS,STA_VEL>{
 public:
  typedef ResidualPosBase<OUT_POS,STA_POS,STA_VEL> Base;
  using Base::dt_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  ResidualPos(): ResidualPosBase<OUT_POS,STA_POS,STA_VEL>(true,true,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_POS>() = cur.template Get<STA_POS>() - pre.template Get<STA_POS>()
                                  - dt_*pre.template Get<STA_VEL>();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_POS)) = -Mat3::Identity();
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_VEL)) = -Mat3::Identity()*dt_;
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_POS),pre.Start(STA_POS)) = Mat3::Identity();
    return 0;
  }
};

template<int OUT_VEL, int STA_VEL>
using ResdidualVelBase = Residual<ElementVector<Element<Vec3,OUT_VEL>>,
                                  ElementVector<Element<Vec3,STA_VEL>>,
                                  ElementVector<Element<Vec3,STA_VEL>>,
                                  MeasA>;

template<int OUT_VEL, int STA_VEL>
class ResidualVel: public ResdidualVelBase<OUT_VEL,STA_VEL>{
 public:
  typedef ResdidualVelBase<OUT_VEL,STA_VEL> Base;
  using Base::dt_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  ResidualVel(): ResdidualVelBase<OUT_VEL,STA_VEL>(true,false,true){}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    out.template Get<OUT_VEL>() = cur.template Get<STA_VEL>() - pre.template Get<STA_VEL>() - dt_*meas_->GetAcc();
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_VEL),pre.Start(STA_VEL)) = -Mat3::Identity();
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    J.block<3,3>(Output::Start(OUT_VEL),pre.Start(STA_VEL)) = Mat3::Identity();
    return 0;
  }
};


int main(int argc, char** argv){
  typedef Filter<ResidualPos<0,0,1>,ResidualVel<0,1>> FilterA;
  FilterA filter;
  FilterA::State state_pre;
  state_pre.SetRandom();
  FilterA::State state_cur;
  state_cur.SetRandom();
  filter.JacTestAll(1e-6,1e-8,state_pre,state_cur);
  std::cout << state_cur.Print() << std::endl;

  std::cout << Print(filter.GetMinMaxTime()) << std::endl;
  TimePoint start = Clock::now();
  filter.Update();
  std::cout << Print(filter.GetMinMaxTime()) << std::endl;
  filter.Update();
  filter.AddMeas<1>(start+fromSec(0.01),std::make_shared<MeasA>(NormalRandomNumberGenerator::Instance().GetVec<3>()));
  std::cout << Print(filter.GetMinMaxTime()) << std::endl;
  filter.AddMeas<1>(start+fromSec(0.03),std::make_shared<MeasA>(NormalRandomNumberGenerator::Instance().GetVec<3>()));
  std::cout << Print(filter.GetMinMaxTime()) << std::endl;
  filter.Update();
  std::cout << filter.GetState().Print() << std::endl;
  for(int i=0;i<10;i++){
    filter.AddMeas<1>(start+fromSec(0.05+0.02*i),std::make_shared<MeasA>(NormalRandomNumberGenerator::Instance().GetVec<3>()));
  }
  filter.Update();
  std::cout << filter.PrintConnectivity() << std::endl;

  return 0;
}
