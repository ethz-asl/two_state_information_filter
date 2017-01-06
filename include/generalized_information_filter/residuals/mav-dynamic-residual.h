#ifndef GIF_MAVDYNAMICRESIDUAL_H_
#define GIF_MAVDYNAMICRESIDUAL_H_

#include <Eigen/Dense>

namespace GIF {

/*! \brief Rotor Speed Measurements
 *         ElementVector that can be used to hold rotor speed measurements.
 */
template<int NumRot>
class RotorSpeedMeasurement : public ElementVector {
 public:
  RotorSpeedMeasurement(const Vec<NumRot>& rot = ElementTraits<Vec<NumRot>>::Identity())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
          new ElementPack<Vec<NumRot>>({"rot"}))),
          rot_(ElementVector::GetValue<Vec<NumRot>>("rot")) {
    rot_ = rot;
  }
  Vec<NumRot>& rot_;
};

/*! \brief Rotor Speed Measurements Residual
 *         Based on an aerodynamic model
 */
template<int NumRot>
class MavDynamicResidual : public BinaryResidual<
ElementPack<Vec<6>>,
ElementPack<Vec3, Vec3, Quat, double>,
ElementPack<Vec3, Vec3>,
ElementPack<Vec<6>>,
RotorSpeedMeasurement<NumRot>> {
 public:
  using mtBinaryResidual = BinaryResidual<
      ElementPack<Vec<6>>,
      ElementPack<Vec3, Vec3, Quat, double>,
      ElementPack<Vec3, Vec3>,
      ElementPack<Vec<6>>,
      RotorSpeedMeasurement<NumRot>>;
  using mtBinaryResidual::meas_;
  using mtBinaryResidual::dt_;
  using mtBinaryResidual::PreDefinition;
  using mtBinaryResidual::CurDefinition;
  using mtBinaryResidual::NoiDefinition;
  using mtBinaryResidual::JacFDPre;
  using mtBinaryResidual::JacFDCur;
  using mtBinaryResidual::JacFDNoi;

  static constexpr int kNumRot_ = NumRot;

  MavDynamicResidual(const std::string& name,
                     const std::string& innName = "dyn",
                     const std::array<std::string,4>& preName = {"MvM", "MwM", "qIM", "m"},
                     const std::array<std::string,2>& curName = {"MvM", "MwM"},
                     const std::string& noiName = "dyn")
      : mtBinaryResidual(name,{innName},preName,curName,{noiName},false,false,false), // TODO
        g_(0,0,-9.81){
    calibrateMass_ = true;
    double planeOffset = -0.016;
    double armLength = 0.215;
    for(int i=0; i<NumRot; i++){
      BrBA_[i](0) = cos((30+i*60)/180*M_PI)*armLength;
      BrBA_[i](1) = sin((30+i*60)/180*M_PI)*armLength;
      BrBA_[i](2) = planeOffset;
      sign_[i] = (i%2 == 0);
    }
    cD_ = 0.018529;
    cM_ = 4.2101e-8;
    cT_ = 5.73311e-6;
    m_ = 1.678;
    I_.setIdentity();
    I_(0,0) = 1.48542e-2;
    I_(1,1) = 1.60995e-2;
    I_(2,2) = 1.34006e-2;
    BrBC_ = Vec3(2.40650e-2, -1.39757e-2, -3.71522e-2);
    SetExtrinsics(Vec3(-0.0275, 0.0158, 0.0377),Quat(1,0,0,0));
  }

  virtual ~MavDynamicResidual() {
  }

  void Eval(Vec<6>& dyn_inn, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                             const double& m_pre,
                             const Vec3& MvM_cur, const Vec3& MwM_cur,
                             const Vec<6>& dyn_noi) const {
    // Compute Bforces
    std::array<Vec<3>,NumRot> Bforces;
    for(int i=0; i<NumRot; i++){
      const Vec<3> MrMAi = qMB_.rotate(Vec3(BrBA_[i]-BrBM_));
      const double T = cT_*meas_->rot_(i)*meas_->rot_(i);
      Bforces[i] = -T*cD_*qMB_.inverseRotate(Vec3(MvM_pre + MwM_pre.cross(MrMAi)));
      Bforces[i](2) = T;
    }

    // Compute Bmoments
    std::array<Vec<3>,NumRot> Bmoments;
    for(int i=0; i<NumRot; i++){
      Bmoments[i].setZero();
      Bmoments[i](2) = (sign_[i]*2-1)*cM_*meas_->rot_(i)*meas_->rot_(i);
    }

    // Compute total force and moment (w.r.t. CoM)
    Vec<3> B_F(0,0,0);
    for(int i=0; i<NumRot; i++){
      B_F += Bforces[i];
    }
    Vec<3> B_M(0,0,0);
    for(int i=0; i<NumRot; i++){
      const Vec<3> BrCAi = BrBA_[i]-BrBC_;
      B_M += Bmoments[i] + BrCAi.cross(Bforces[i]);
    }

    // Evaluate residual by means of equation of motion
    const Vec3 MrMC = qMB_.rotate(Vec3(BrBC_ - BrBM_));
    const Vec3 dMvM = (MvM_cur - MvM_pre)/dt_;
    const Vec3 dMwM = (MwM_cur - MwM_pre)/dt_;
    dyn_inn.head<3>() = 1/m_*qMB_.rotate(B_F) + qIM_pre.inverseRotate(g_) - MwM_pre.cross(MvM_pre) - MwM_pre.cross(MwM_pre.cross(MrMC)) - dMwM.cross(MrMC) - dMvM + dyn_noi.head<3>()*sqrt(dt_);
    const Vec3 BwM_pre = qMB_.inverseRotate(MwM_pre);
    const Vec3 dBwM = qMB_.inverseRotate(dMwM);
    dyn_inn.tail<3>() = B_M - BwM_pre.cross(I_*BwM_pre) - I_*dBwM + dyn_noi.tail<3>()*sqrt(dt_);
  }
  void JacPre(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<6>& dyn_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<Vec<6>>(0) = dyn_noi;
    JacFDPre(J,pre,cur,noi,1e-6);
  }
  void JacCur(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<6>& dyn_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<Vec<6>>(0) = dyn_noi;
    JacFDCur(J,pre,cur,noi,1e-6);
  }
  void JacNoi(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<6>& dyn_noi) const {
    J.setIdentity();
    J = J*sqrt(dt_);
  }
  void SetExtrinsics(Vec3 BrBM, Quat qMB){
    BrBM_  = BrBM;
    qMB_ = qMB;
  }

 protected:
  Vec3 BrBM_;
  Vec3 BrBC_; // Center of mass
  Quat qMB_;
  std::array<Vec<3>,NumRot> BrBA_;
  std::array<bool,NumRot> sign_;
  double cD_;
  double cM_;
  double cT_;
  double m_;
  Vec3 g_;
  Mat3 I_;
  bool calibrateMass_;
};


}

#endif /* GIF_MAVDYNAMICRESIDUAL_H_ */
