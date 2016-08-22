/*
 * IMUPrediction.hpp
 *
 *  Created on: Aug 22, 2016
 *      Author: Bloeschm
 */

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Prediction.hpp"

namespace GIF{

class IMUMeas: public State{
 public:
  IMUMeas(const V3D& gyr = V3D(0,0,0), const V3D& acc = V3D(0,0,0)):
      State(ElementPack<V3D,V3D>::makeStateDefinition({"gyr","acc"})),
      gyr_(State::getValue<V3D>("gyr")),
      acc_(State::getValue<V3D>("acc")){
    gyr_ = gyr;
    acc_ = acc;
  };
  V3D& gyr_;
  V3D& acc_;
};

class IMUPrediction: public Prediction<ElementPack<V3D,V3D,V3D,V3D,QPD>,ElementPack<V3D,V3D,V3D,V3D,V3D>,IMUMeas>{
 public:
  IMUPrediction(): mtPrediction({"pos","vel","gyb","acb","att"},{"pos","vel","gyb","acb","att"}), g_(0,0,-9.81){
    dt_ = 0.1;
    meas_.reset(new IMUMeas());
  };
  virtual ~IMUPrediction(){};
  void evalPredictionImpl(        V3D& posPos ,      V3D& velPos  ,      V3D& gybPos  ,      V3D& acbPos  ,      QPD& attPos,
                            const V3D& posPre ,const V3D& velPre  ,const V3D& gybPre  ,const V3D& acbPre  ,const QPD& attPre,
                            const V3D& posNoi ,const V3D& velNoi  ,const V3D& gybNoi  ,const V3D& acbNoi  ,const V3D& attNoi) const{
    const V3D gyr = meas_->gyr_ - gybPre + attNoi/sqrt(dt_);
    const V3D acc = meas_->acc_ - acbPre + velNoi/sqrt(dt_);
    const V3D dOmega = dt_*gyr;
    QPD dQ = dQ.exponentialMap(dOmega);
    posPos = posPre + dt_*(attPre.rotate(velPre) + posNoi/sqrt(dt_));
    velPos = (M3D::Identity() - gSM(dOmega))*velPre + dt_*(acc-attPre.inverseRotate(g_));
    gybPos = gybPre + gybNoi*sqrt(dt_);
    acbPos = acbPre + acbNoi*sqrt(dt_);
    attPos = attPre*dQ;
  }
  void jacPrePredictionImpl(  MXD& J,
                              const V3D& posPre ,const V3D& velPre  ,const V3D& gybPre  ,const V3D& acbPre  ,const QPD& attPre,
                              const V3D& posNoi ,const V3D& velNoi  ,const V3D& gybNoi  ,const V3D& acbNoi  ,const V3D& attNoi) const{
    J.setZero();
    const V3D gyr = meas_->gyr_ - gybPre + attNoi/sqrt(dt_);
    const V3D acc = meas_->acc_ - acbPre + velNoi/sqrt(dt_);
    const V3D dOmega = dt_*gyr;
    setJacBlockPre<POS,POS>(J,M3D::Identity());
    setJacBlockPre<POS,VEL>(J,dt_*MPD(attPre).matrix());
    setJacBlockPre<POS,ATT>(J,dt_*gSM(attPre.rotate(velPre)));
    setJacBlockPre<VEL,VEL>(J,(M3D::Identity() - gSM(dOmega)));
    setJacBlockPre<VEL,GYB>(J,-dt_*gSM(gyr));
    setJacBlockPre<VEL,ACB>(J,-dt_*M3D::Identity());
    setJacBlockPre<VEL,ATT>(J,-dt_*MPD(attPre).matrix().transpose()*gSM(g_));
    setJacBlockPre<GYB,GYB>(J,M3D::Identity());
    setJacBlockPre<ACB,ACB>(J,M3D::Identity());
    setJacBlockPre<ATT,GYB>(J,-dt_*MPD(attPre).matrix()*Lmat(dOmega));
    setJacBlockPre<ATT,ATT>(J,M3D::Identity());
  }
  void jacNoiPredictionImpl(  MXD& J,
                              const V3D& posPre ,const V3D& velPre  ,const V3D& acbPre  ,const V3D& gybPre  ,const QPD& attPre,
                              const V3D& posNoi ,const V3D& velNoi  ,const V3D& acbNoi  ,const V3D& gybNoi  ,const V3D& attNoi) const{
    J.setZero();
    const V3D gyr = meas_->gyr_ - gybPre + attNoi/sqrt(dt_);
    const V3D acc = meas_->acc_ - acbPre + velNoi/sqrt(dt_);
    const V3D dOmega = dt_*gyr;
    setJacBlockNoi<POS,POS>(J,sqrt(dt_)*M3D::Identity());
    setJacBlockNoi<VEL,VEL>(J,sqrt(dt_)*M3D::Identity());
    setJacBlockNoi<VEL,ATT>(J,sqrt(dt_)*gSM(velPre));
    setJacBlockNoi<GYB,GYB>(J,sqrt(dt_)*M3D::Identity());
    setJacBlockNoi<ACB,ACB>(J,sqrt(dt_)*M3D::Identity());
    setJacBlockNoi<ATT,ATT>(J,sqrt(dt_)*MPD(attPre).matrix()*Lmat(dOmega));
  }

 protected:
  double dt_;
  const V3D g_;
  enum Elements{POS, VEL, GYB, ACB, ATT};
};

}
