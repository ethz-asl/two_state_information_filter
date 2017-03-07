#ifndef GIF_LEGDYNAMICRESIDUAL_H_
#define GIF_LEGDYNAMICRESIDUAL_H_

#include <Eigen/Dense>

namespace GIF {

/*! \brief Leg dynamics Measurement (torques)
 *         ElementVector that can be used to hold leg dynamics measurements. It is composed of an
 *         array of vectors containing the single torque measurements.
 */
template<int NumLeg, int NumDof>
class DynamicMeasurement : public ElementVector {
 public:
  typedef std::array<Vec<NumDof>,NumLeg> Arr;
  DynamicMeasurement(const Arr& enc = ElementTraits<Arr>::Identity(),
                     const Arr& end = ElementTraits<Arr>::Identity(),
                     const Arr& dyn = ElementTraits<Arr>::Identity())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
          new ElementPack<Arr,Arr,Arr>({"enc", "end", "dyn"}))),
          enc_(ElementVector::GetValue<Arr>("enc")),
          end_(ElementVector::GetValue<Arr>("end")),
          dyn_(ElementVector::GetValue<Arr>("dyn")) {
    enc_ = enc;
    end_ = end;
    dyn_ = dyn;
    for(int i=0;i<NumLeg;i++){
      contact_flag_[i] = false;
    }
  }
  Arr& enc_;
  Arr& end_;
  Arr& dyn_;
  std::array<bool,NumLeg> contact_flag_;
};

/*! \brief Leg Dynamics Binary Residual
 *         Integrates information derived from the equations of motions together with joint torque
 *         measurements.
 */
template<typename RobotModel>
class LegDynamicResidual : public BinaryResidual<
    ElementPack<Vec<RobotModel::kNumLeg*RobotModel::kNumDof+6>>,
    ElementPack<Vec3, Vec3, Quat, double, Vec3>,
    ElementPack<Vec3, Vec3>,
    ElementPack<Vec<RobotModel::kNumLeg*RobotModel::kNumDof+6>>,
    DynamicMeasurement<RobotModel::kNumLeg, RobotModel::kNumDof>> {
 public:
  static constexpr int kDof_ = RobotModel::kNumLeg*RobotModel::kNumDof+6;
  typedef Vec<kDof_> DynRes;
  using mtBinaryResidual = BinaryResidual<
      ElementPack<DynRes>,
      ElementPack<Vec3, Vec3, Quat, double, Vec3>,
      ElementPack<Vec3, Vec3>,
      ElementPack<Vec<kDof_>>,
      DynamicMeasurement<RobotModel::kNumLeg, RobotModel::kNumDof>>;
  using mtBinaryResidual::meas_;
  using mtBinaryResidual::dt_;
  using mtBinaryResidual::PreDefinition;
  using mtBinaryResidual::CurDefinition;
  using mtBinaryResidual::NoiDefinition;
  using mtBinaryResidual::JacFDPre;
  using mtBinaryResidual::JacFDCur;
  using mtBinaryResidual::JacFDNoi;
  enum Elements {DYN};

  LegDynamicResidual(const std::string& name,
                     const std::string& innName = "dyn",
                     const std::array<std::string,5>& preName = {"MvM", "MwM", "qIM", "m", "mo"},
                     const std::array<std::string,2>& curName = {"MvM", "MwM"},
                     const std::string& noiName = "dyn")
      : mtBinaryResidual(name,{innName},preName,curName,{noiName},false,false,false) { // TODO
    hasPreviousEncData_ = false;
    verbose_ = true;
    calibrateMass_ = true;
    calibrateMassOffset_ = true;
  }

  virtual ~LegDynamicResidual() {
  }

  void Eval(DynRes& dyn_inn, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                             const double& m_pre, const Vec3& mo_pre,
                             const Vec3& MvM_cur, const Vec3& MwM_cur,
                             const DynRes& dyn_noi) const {
    if(hasPreviousEncData_){
      // Change to body coordinates
      const Vec3 BvB_pre = qMB_.inverseRotate(Vec3(MvM_pre + MwM_pre.cross(qMB_.rotate(BrBM_))));
      const Vec3 BwB_pre = qMB_.inverseRotate(MwM_pre);
      const Quat qIB_pre = qIM_pre*qMB_;
      const Vec3 BvB_cur = qMB_.inverseRotate(Vec3(MvM_cur + MwM_cur.cross(qMB_.rotate(BrBM_))));
      const Vec3 BwB_cur = qMB_.inverseRotate(MwM_cur);

      // Get generalized quantities
      genPos_pre = model_->computeGeneralizedPositions(Vec3(0,0,0),qIB_pre,enc_pre_);
      genVel_pre = model_->computeGeneralizedVelocities(qIB_pre.rotate(BvB_pre),
                                                        BwB_pre,end_pre_);
      Quat dQ = dQ.exponentialMap(dt_ * BwB_cur);
      genVel_cur = model_->computeGeneralizedVelocities((qIB_pre*dQ).rotate(BvB_cur),
                                                        BwB_cur,meas_->end_);
      genAcc = (genVel_cur-genVel_pre)/dt_;

      // Pass to model
      if(calibrateMass_){
        model_->setMass(m_pre);
      }
      if(calibrateMassOffset_){
        model_->setMassOffset(mo_pre);
      }
      model_->setGeneralizedPositions(Vec3(0,0,0),qIB_pre,enc_pre_);
      model_->setGeneralizedVelocities(qIB_pre.rotate(BvB_pre),BwB_pre,end_pre_);
      model_->updateKinematics(true,true,false);

      // Get dynamic equation quantities
      MatX M((int)kDof_,(int)kDof_);
      M.setZero();
      model_->getMassInertiaMatrix(M);
      VecX h((int)kDof_);
      h.setZero();
      model_->getNonlinearEffects(h);
      MatX S((int)kDof_-6,(int)kDof_);
      S.setZero();
      model_->getSelectionMatrix(S);
      MatX J(RobotModel::kNumLeg*3,(int)kDof_);
      for(int i=0;i<RobotModel::kNumLeg;i++){
        MatX J_sub = model_->getJacobianTranslationWorldToFoot(i);
        J.block<3,kDof_>(3*i,0) = J_sub;
      }
      VecX jointTorques((int)kDof_-6);
      for(int i=0;i<RobotModel::kNumLeg;i++){
        jointTorques.segment<RobotModel::kNumDof>(RobotModel::kNumDof*i) = meas_->dyn_[i];
      }

      // Make projection
      /*
      J^T = Q R P^-1 = [Q1, Q2] [R1^T, 0^T]^T  P^-1 = Q1 R1 P^-1
      N = Q2^T
      N*J^T = Q2^T Q1 R1 P^-1 = 0
      */
      Eigen::ColPivHouseholderQR<MatX> qr(J.transpose());
      int rank = qr.rank();
      MatX N = MatX(qr.householderQ()).block(0,rank,18,18-rank).transpose();

      // Compute residual
      VecX y(18);
      y.setZero();
      y.head(18-rank) = N*(M*genAcc+h-S.transpose()*jointTorques);
//      std::cout << y.head(18-rank).transpose() << std::endl;

//      if(verbose_){
//        std::cout << "----" << std::endl;
//        std::cout << qIB_pre << std::endl;
//        std::cout << (enc_pre_[0]).transpose() << std::endl;
//        std::cout << (enc_pre_[1]).transpose() << std::endl;
//        std::cout << (enc_pre_[2]).transpose() << std::endl;
//        std::cout << (enc_pre_[3]).transpose() << std::endl;
//        std::cout << genPos_pre.transpose() << std::endl;
//        std::cout << (qIB_pre.rotate(BvB_pre)).transpose() << std::endl;
//        std::cout << (BwB_pre).transpose() << std::endl;
//        std::cout << (end_pre_[0]).transpose() << std::endl;
//        std::cout << (end_pre_[1]).transpose() << std::endl;
//        std::cout << (end_pre_[2]).transpose() << std::endl;
//        std::cout << (end_pre_[3]).transpose() << std::endl;
//        std::cout << genVel_pre.transpose() << std::endl;
//        std::cout << genVel_cur.transpose() << std::endl;
//        std::cout << genAcc.transpose() << std::endl;
//        std::cout << "----" << std::endl;
//        std::cout << "J: " << std::endl << J << std::endl;
//        std::cout << "----" << std::endl;
//        std::cout << jointTorques.transpose() << std::endl;
//        std::cout << "----" << std::endl;
//        std::cout << "Ma: " << (N*M*genAcc).transpose() << std::endl;
//        std::cout << "h : " << (N*h).transpose() << std::endl;
//        std::cout << "St: " << (-N*S.transpose()*jointTorques).transpose() << std::endl;
//        std::cout << "T : " << y.transpose() << std::endl;
//
//
//
//        Vec3 pos1;
//        model_->getPositionWorldToBody(pos1,model_->getBranchEnumFromLimbUInt(0),
//                                       quadruped_model::BodyNodeEnum::FOOT);
//        std::cout << "Foothold1: " << pos1.transpose() << std::endl;
//        Vec3 dif(-0.3,0.2,0.1);
//        double d =  1e-8;
//        std::array<Vec<RobotModel::kNumDof>,RobotModel::kNumLeg> enc_pre;
//        enc_pre[0] = enc_pre_[0]+d*dif;
//        enc_pre[1] = enc_pre_[1];
//        enc_pre[2] = enc_pre_[2];
//        enc_pre[3] = enc_pre_[3];
//        model_->setGeneralizedPositions(Vec3(0,0,0),qIB_pre,enc_pre);
//        model_->updateKinematics(true,true,false);
//        Vec3 pos2;
//        model_->getPositionWorldToBody(pos2,model_->getBranchEnumFromLimbUInt(0),
//                                       quadruped_model::BodyNodeEnum::FOOT);
//        std::cout << "Foothold2: " << pos2.transpose() << std::endl;
//        std::cout << "Gradient1: " << (pos2-pos1).transpose()/d << std::endl;
//        std::cout << "Gradient2: " << (J.block<3,3>(0,6)*dif).transpose() << std::endl;
//
//
//        Quat dQQ = dQQ.exponentialMap(d * dif);
//        Quat qIB_pre_perturbed = qIB_pre*dQQ;
//        model_->setGeneralizedPositions(Vec3(0,0,0),qIB_pre_perturbed,enc_pre_);
//        model_->updateKinematics(true,true,false);
//        model_->getPositionWorldToBody(pos2,model_->getBranchEnumFromLimbUInt(0),
//                                       quadruped_model::BodyNodeEnum::FOOT);
//        std::cout << "Foothold3: " << pos2.transpose() << std::endl;
//        std::cout << "Gradient3: " << (pos2-pos1).transpose()/d << std::endl;
//        std::cout << "Gradient4: " << (J.block<3,3>(0,3)*dif).transpose() << std::endl;
//
//
//
//      }
      dyn_inn = y + dyn_noi*sqrt(dt_);
    } else {
      dyn_inn = dyn_noi*sqrt(dt_);
    }
  }
  void JacPre(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const Vec3& mo_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const DynRes& dyn_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    pre.GetValue<Vec3>(4) = mo_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<DynRes>(0) = dyn_noi;
    verbose_ = false;
    JacFDPre(J,pre,cur,noi,1e-6);
    J.col(8).setZero(); // Ensure that Jacobian w.r.t. equals zero
    verbose_ = true;
  }
  void JacCur(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const Vec3& mo_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const DynRes& dyn_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    pre.GetValue<Vec3>(4) = mo_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<DynRes>(0) = dyn_noi;
    verbose_ = false;
    JacFDCur(J,pre,cur,noi,1e-6);
    verbose_ = true;
  }
  void JacNoi(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const Vec3& mo_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const DynRes& dyn_noi) const {
    J.setIdentity();
    J = J*sqrt(dt_);
  }

  void SetModelPtr(const std::shared_ptr<RobotModel>& model) {
    model_ = model;
  }
  void SetPreviousEncoderData(){
    enc_pre_ = meas_->enc_;
    end_pre_ = meas_->end_;
    hasPreviousEncData_ = true;
  }
  void ResetPreviousEncoderData(){
    hasPreviousEncData_ = false;
    for(int i=0; i<RobotModel::kNumLeg; i++){
      enc_pre_[i].setZero();
      end_pre_[i].setZero();
    }
  }

  void SetExtrinsics(Vec3 BrBM, Quat qMB){
    BrBM_  = BrBM;
    qMB_ = qMB;
  }
  mutable bool verbose_;

 protected:
  std::shared_ptr<RobotModel> model_;
  Vec3 BrBM_;
  Quat qMB_;
  std::array<Vec<RobotModel::kNumDof>,RobotModel::kNumLeg> enc_pre_;
  std::array<Vec<RobotModel::kNumDof>,RobotModel::kNumLeg> end_pre_;
  bool hasPreviousEncData_;
  mutable VecX genPos_pre;
  mutable VecX genVel_pre;
  mutable VecX genVel_cur;
  mutable VecX genAcc;
  bool calibrateMass_;
  bool calibrateMassOffset_;
};


}

#endif /* GIF_LEGDYNAMICRESIDUAL_H_ */
