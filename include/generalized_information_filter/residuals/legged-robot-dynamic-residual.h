#ifndef GIF_LEGGEDROBOTDYNAMICRESIDUAL_HPP_
#define GIF_LEGGEDROBOTDYNAMICRESIDUAL_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"

namespace GIF {

/*! \brief Leg Kinematic and Dynamic Measurement
 *         ElementVector that can be used to hold leg kinematic and dynamic measurements. It is
 *         composed of three arrays of vectors containing the single measurements
 */
template<int NumLeg, int NumDof>
class DynamicMeasurement : public ElementVector {
 public:
  typedef std::array<Vec<NumDof>,NumLeg> Kin;
  DynamicMeasurement(const Kin& ang = ElementTraits<Kin>::Identity(),
                     const Kin& vel = ElementTraits<Kin>::Identity(),
                     const Kin& tor = ElementTraits<Kin>::Identity())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
          new ElementPack<Kin,Kin,Kin>({"ang","vel","tor"}))),
        ang_(ElementVector::GetValue<Kin>("ang")),
        vel_(ElementVector::GetValue<Kin>("vel")),
        tor_(ElementVector::GetValue<Kin>("tor")) {
    ang_ = ang;
    vel_ = vel;
    tor_ = tor;
    for(int i=0;i<NumLeg;i++){
      contact_flag_[i] = false;
    }
  }
  Kin& ang_;
  Kin& vel_;
  Kin& tor_;
  std::array<bool,NumLeg> contact_flag_;
};

/*! \brief Legged robot dynamic prediction
 *         Leverages the dynamic measurements (joint torques) in order to related to subsequent
 *         state. It employes the projected equation of motion for deriving error terms on the
 *         acceleration level.
 *
 *         Coordinate frames:
 *           B: Body
 *           M: IMU
 *           F: Foothold
 */
template<typename KinModel>
class LeggedRobotDynamicResidual : public BinaryResidual<
                                        ElementPack<Vec<KinModel::kNumLeg*KinModel::kNumDof+6>>,
                                        ElementPack<Vec3,Vec3,Vec3,Quat>,
                                        ElementPack<Vec3,Vec3,Vec3,Quat>,
                                        ElementPack<Vec<KinModel::kNumLeg*KinModel::kNumDof+6>>,
                                        DynamicMeasurement<KinModel::kNumLeg,
                                                           KinModel::kNumDof>> {
 public:
  static constexpr int kNumDof = KinModel::kNumLeg*KinModel::kNumDof+6;
  typedef std::array<std::string,4> String4;
  typedef std::array<std::string,1> String1;
  using mtBinaryResidual = BinaryResidual<ElementPack<Vec<kNumDof>>,
                                          ElementPack<Vec3,Vec3,Vec3,Quat>,
                                          ElementPack<Vec3,Vec3,Vec3,Quat>,
                                          ElementPack<Vec<kNumDof>>,
                                          DynamicMeasurement<KinModel::kNumLeg,
                                                             KinModel::kNumDof>>;
  using mtBinaryResidual::meas_;
  enum Elements {POS,VEL,ROR,ATT};

  LeggedRobotDynamicResidual(const String1& innName = {"dyn"},
                               const String4& preName = {"IrIM","MvM","MwM","qIM"},
                               const String4& curName = {"IrIM","MvM","MwM","qIM"},
                               const String1& noiName = {"dyn"})
      : mtBinaryResidual(innName, preName, curName, noiName, false, true, true) {
  }

  virtual ~LeggedRobotDynamicResidual() {
  }

  void Eval(Vec<kNumDof>& inn,
            const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
            const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur, const Quat& qIM_cur,
            const Vec<kNumDof>& noi) const {
//    typename DynamicModelMinimal<nFeet,estFP,KinModel>::mtInnovation dynResidual;
//    typename DynamicModelMinimal<nFeet,estFP,KinModel>::mtNoise dynNoise;
//    dynNoise.template get<0>() = noise.template get<mtNoise::_dyn>();
//    dynamicModelMinimal_.evalResidual(dynResidual,previousState,currentState,dynNoise,dt);
//    y.template get<mtInnovation::_dyn>() = dynResidual.template get<0>();

    // Compute generalized positions and velocities for previous and current state


//    mpCommonModel_->getGeneralizedPositions(genPosPrev_,previousState);
//    mpCommonModel_->getGeneralizedPositions(genPosCurr_,currentState);
//    mpCommonModel_->getGeneralizedVelocities(genVelPrev_,previousState);
//    mpCommonModel_->getGeneralizedVelocities(genVelCurr_,currentState);
//    genAcc_ = (genVelCurr_ - genVelPrev_)/dt;
//    mpCommonModel_->getMassInertiaMatrix(M_,genPosPrev_);
//    mpCommonModel_->getNonlinearEffects(h_,genPosPrev_,genVelPrev_);
//    mpCommonModel_->getFullContactJacobian(J_,genPosPrev_);
//    mpCommonModel_->getSelectionMatrix(S_);
//    previousState.aux().meas_.getJointTorques(jointTorques_);
//    mpCommonModel_->getContactForcesInWorld(contactForces_,previousState);
//    N_ = Eigen::Matrix<double,nFeet*3+6,nFeet*3+6>::Identity() - J_.transpose()*(J_*J_.transpose()).inverse()*J_;
////    std::cout << "=======================================" << std::endl;
//    y.template get<0>() = N_*(M_*genAcc_+h_-S_.transpose()*jointTorques_) + noise.template get<0>()/sqrt(dt); // TODO: check noise, evaluate magnitude, evaluate different residuals
  }
  void JacPre(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur, const Quat& qIM_cur,
              const Vec<kNumDof>& noi) const {
    J.setZero();
  }
  void JacCur(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur, const Quat& qIM_cur,
              const Vec<kNumDof>& noi) const {
    J.setZero();
  }
  void JacNoi(MatX& J,
              const Vec3& IrIM_pre, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
              const Vec3& IrIM_cur, const Vec3& MvM_cur, const Vec3& MwM_cur, const Quat& qIM_cur,
              const Vec<kNumDof>& noi) const {
    J.setZero();
  }

  void SetModelPtr(const std::shared_ptr<KinModel>& model) {
    model_ = model;
  }

  void SetExtrinsics(Vec3 BrBM, Quat qMB){
    BrBM_  = BrBM;
    qMB_ = qMB;
  }

 protected:
  std::shared_ptr<KinModel> model_;
  Vec3 BrBM_;
  Quat qMB_;
};

}

#endif /* GIF_LEGGEDROBOTDYNAMICRESIDUAL_HPP_ */
