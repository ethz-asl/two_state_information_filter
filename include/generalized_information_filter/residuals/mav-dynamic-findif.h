#ifndef GIF_MAVDYNAMICFINDIF_H_
#define GIF_MAVDYNAMICFINDIF_H_

#include <Eigen/Dense>
#include "generalized_information_filter/binary-residual.h"

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
class MavDynamicFindif : public BinaryResidual<
    ElementPack<Vec<3>,Vec<3>>,
    ElementPack<Vec3, Vec3, Quat, double, double, double, double, Vec3, Vec3, Vec3, Quat>,
    ElementPack<Vec3, Vec3>,
    ElementPack<Vec<3>,Vec<3>>,
RotorSpeedMeasurement<NumRot>> {
 public:
  using mtBinaryResidual = BinaryResidual<
      ElementPack<Vec<3>,Vec<3>>,
      ElementPack<Vec3, Vec3, Quat, double, double, double, double, Vec3, Vec3, Vec3, Quat>,
      ElementPack<Vec3, Vec3>,
      ElementPack<Vec<3>,Vec<3>>,
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

  MavDynamicFindif(const std::string& name,
                   const std::array<std::string,2>& innName = {"dyn_pos", "dyn_att"},
                   const std::array<std::string,11>& preName = {"MvM", "MwM", "qIM", "m", "cT",
                                                                "cM", "cD", "BrBM", "BrBC", "I", "qMB"},
                   const std::array<std::string,2>& curName = {"MvM", "MwM"},
                   const std::array<std::string,2>& noiName = {"dyn_pos", "dyn_att"})
      : mtBinaryResidual(name,innName,preName,curName,noiName,false,true,true),
        g_(0,0,-9.81){
    SetCalibrationFlags(0,0,0,0,0,0,0,0);
    double planeOffset = -0.016;
    double armLength = 0.215;
    for(int i=0; i<NumRot; i++){
      BrBA_[i](0) = cos((30+i*60)/180*M_PI)*armLength;
      BrBA_[i](1) = sin((30+i*60)/180*M_PI)*armLength;
      BrBA_[i](2) = planeOffset;
      sign_[i] = (i%2 == 0);
    }
    SetMass(1.678);
    SetAerodynamicCoefficient(5.73311e-6*1e6,-4.2101e-8*1e6,0.018529);
    SetInertiaDiagonal(1.48542e-2,1.60995e-2,1.34006e-2);
    SetComOffset(Vec3(0.0241, -0.0140, -0.0372));
    SetImuOffset(Vec3(-0.0275, 0.0158, 0.0377),Quat(1,0,0,0));
  }

  virtual ~MavDynamicFindif() {
  }

  void Eval(Vec<3>& dyn_pos_inn, Vec<3>& dyn_att_inn,
            const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
            const double& m_pre, const double& cT_pre, const double& cM_pre, const double& cD_pre,
            const Vec3& BrBM_pre, const Vec3& BrBC_pre, const Vec3& I_pre, const Quat& qMB_pre,
            const Vec3& MvM_cur, const Vec3& MwM_cur,
            const Vec<3>& dyn_pos_noi, const Vec<3>& dyn_att_noi) const {
    double m = calibrationFlags_[0] ? m_pre : m_;
    double cT = calibrationFlags_[1] ? cT_pre : cT_;
    double cM = calibrationFlags_[2] ? cM_pre : cM_;
    double cD = calibrationFlags_[3] ? cD_pre : cD_;
    Vec3 BrBM = calibrationFlags_[4] ? BrBM_pre : BrBM_;
    Vec3 BrBC = calibrationFlags_[5] ? BrBC_pre : BrBC_;
    Mat3 I = calibrationFlags_[6] ? I_pre.asDiagonal() : I_;
    Quat qMB = calibrationFlags_[7] ? qMB_pre : qMB_;

    // Compute Bforces
    std::array<Vec<3>,NumRot> Bforces;
    for(int i=0; i<NumRot; i++){
      const Vec<3> MrMAi = qMB.rotate(Vec3(BrBA_[i]-BrBM));
      const double T = cT*1e-6*meas_->rot_(i)*meas_->rot_(i);
      Bforces[i] = -T*cD*qMB.inverseRotate(Vec3(MvM_pre + MwM_pre.cross(MrMAi)));
      Bforces[i](2) = T;
    }

    // Compute Bmoments
    std::array<Vec<3>,NumRot> Bmoments;
    for(int i=0; i<NumRot; i++){
      Bmoments[i].setZero();
      Bmoments[i](2) = (sign_[i]*2-1)*cM*1e-6*meas_->rot_(i)*meas_->rot_(i);
    }

    // Compute total force and moment (w.r.t. CoM)
    Vec<3> B_F(0,0,0);
    for(int i=0; i<NumRot; i++){
      B_F += Bforces[i];
    }
    Vec<3> B_M(0,0,0);
    for(int i=0; i<NumRot; i++){
      const Vec<3> BrCAi = BrBA_[i]-BrBC;
      B_M += Bmoments[i] + BrCAi.cross(Bforces[i]);
    }

    // Evaluate residual by means of equation of motion
    const Vec3 MrMC = qMB.rotate(Vec3(BrBC - BrBM));
    const Vec3 dMvM = (MvM_cur - MvM_pre)/dt_;
    const Vec3 dMwM = (MwM_cur - MwM_pre)/dt_;
    dyn_pos_inn = (1/m*qMB.rotate(B_F) + qIM_pre.inverseRotate(g_) - MwM_pre.cross(MvM_pre)
        - MwM_pre.cross(MwM_pre.cross(MrMC)) - dMwM.cross(MrMC) - dMvM)*dt_ + dyn_pos_noi*sqrt(dt_);
    const Vec3 BwM_pre = qMB.inverseRotate(MwM_pre);
    const Vec3 dBwM = qMB.inverseRotate(dMwM);
    dyn_att_inn = (B_M - BwM_pre.cross(I*BwM_pre) - I*dBwM)*dt_ + dyn_att_noi*sqrt(dt_);
  }
  void JacPre(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const double& cT_pre, const double& cM_pre, const double& cD_pre,
                       const Vec3& BrBM_pre, const Vec3& BrBC_pre, const Vec3& I_pre, const Quat& qMB_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<3>& dyn_pos_noi, const Vec<3>& dyn_att_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    pre.GetValue<double>(4) = cT_pre;
    pre.GetValue<double>(5) = cM_pre;
    pre.GetValue<double>(6) = cD_pre;
    pre.GetValue<Vec3>(7) = BrBM_pre;
    pre.GetValue<Vec3>(8) = BrBC_pre;
    pre.GetValue<Vec3>(9) = I_pre;
    pre.GetValue<Quat>(10) = qMB_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<Vec<3>>(0) = dyn_pos_noi;
    noi.GetValue<Vec<3>>(1) = dyn_att_noi;
    JacFDPre(J,pre,cur,noi,1e-6);
  }
  void JacCur(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const double& cT_pre, const double& cM_pre, const double& cD_pre,
                       const Vec3& BrBM_pre, const Vec3& BrBC_pre, const Vec3& I_pre, const Quat& qMB_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<3>& dyn_pos_noi, const Vec<3>& dyn_att_noi) const {
    ElementVector pre(PreDefinition());
    pre.GetValue<Vec3>(0) = MvM_pre;
    pre.GetValue<Vec3>(1) = MwM_pre;
    pre.GetValue<Quat>(2) = qIM_pre;
    pre.GetValue<double>(3) = m_pre;
    pre.GetValue<double>(4) = cT_pre;
    pre.GetValue<double>(5) = cM_pre;
    pre.GetValue<double>(6) = cD_pre;
    pre.GetValue<Vec3>(7) = BrBM_pre;
    pre.GetValue<Vec3>(8) = BrBC_pre;
    pre.GetValue<Vec3>(9) = I_pre;
    pre.GetValue<Quat>(10) = qMB_pre;
    ElementVector cur(CurDefinition());
    cur.GetValue<Vec3>(0) = MvM_cur;
    cur.GetValue<Vec3>(1) = MwM_cur;
    ElementVector noi(NoiDefinition());
    noi.GetValue<Vec<3>>(0) = dyn_pos_noi;
    noi.GetValue<Vec<3>>(1) = dyn_att_noi;
    JacFDCur(J,pre,cur,noi,1e-6);
  }
  void JacNoi(MatX& J, const Vec3& MvM_pre, const Vec3& MwM_pre, const Quat& qIM_pre,
                       const double& m_pre, const double& cT_pre, const double& cM_pre, const double& cD_pre,
                       const Vec3& BrBM_pre, const Vec3& BrBC_pre, const Vec3& I_pre, const Quat& qMB_pre,
                       const Vec3& MvM_cur, const Vec3& MwM_cur,
                       const Vec<3>& dyn_pos_noi, const Vec<3>& dyn_att_noi) const {
    J.setIdentity();
    J = J*sqrt(dt_);
  }
  void SetImuOffset(Vec3 BrBM, Quat qMB){ // TODO: flip MB
    BrBM_  = BrBM;
    qMB_ = qMB;
  }
  void SetComOffset(Vec3 BrBC){
    BrBC_ = BrBC;
  }
  void SetMass(double m){
    m_ = m;
  }
  void SetAerodynamicCoefficient(double cT, double cM, double cD){
    cT_ = cT;
    cM_ = cM;
    cD_ = cD;
  }
  void SetInertiaDiagonal(double Ix, double Iy, double Iz){
    I_.setIdentity();
    I_(0,0) = Ix;
    I_(1,1) = Iy;
    I_(2,2) = Iz;
  }
  void SetCalibrationFlags(bool m, bool cT, bool cM, bool cD, bool BrBM, bool BrBC, bool I, bool qMB){
    calibrationFlags_[0] = m;
    calibrationFlags_[1] = cT;
    calibrationFlags_[2] = cM;
    calibrationFlags_[3] = cD;
    calibrationFlags_[4] = BrBM;
    calibrationFlags_[5] = BrBC;
    calibrationFlags_[6] = I;
    calibrationFlags_[7] = qMB;
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
  std::array<bool,8> calibrationFlags_;
};


}

#endif /* GIF_MAVDYNAMICFINDIF_H_ */
