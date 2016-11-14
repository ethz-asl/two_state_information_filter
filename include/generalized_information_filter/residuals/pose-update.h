#ifndef GIF_POSEUPDATE_HPP_
#define GIF_POSEUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"

namespace GIF {

/*! \brief Pose Measurement
 *         ElementVector that can be used to hold a generic pose measurements (position + attitude).
 */
class PoseMeas : public ElementVector {
 public:
  PoseMeas(const Vec3& JrJC = Vec3(0, 0, 0), const Quat& qJC = Quat())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Quat>({ "JrJC", "qJC" }))),
        JrJC_(ElementVector::GetValue<Vec3>("JrJC")),
        qJC_(ElementVector::GetValue<Quat>("qJC")) {
    JrJC_ = JrJC;
    qJC_ = qJC;
  }
  Vec3& JrJC_;
  Quat& qJC_;
};

/*! \brief Pose Update
 *         Builds a residual between current estimate pose and an external measured pose. The state
 *         is parametrized by the position (IrIB) and attitude (qIB). Additionally, the offset
 *         between the inertial frame is coestimate (IrIJ and qIJ). The extrinsics calibration of
 *         the pose measurement are assumed to be known  (BrBC and qBC).
 *
 *         Coordinate frames:
 *           B: Body
 *           C: Camera
 *           J: Pose measurement reference inertial frame
 *           I: Estimation reference inertial frame
 */
class PoseUpdate : public UnaryUpdate<ElementPack<Vec3, Quat>,
    ElementPack<Vec3, Quat, Vec3, Quat>, ElementPack<Vec3, Vec3>, PoseMeas> {
 public:
  PoseUpdate(const std::string& name,
             const std::array<std::string,2>& errorName = {"JrJC", "qJC"},
             const std::array<std::string,4>& stateName = {"IrIB", "qIB", "IrIJ", "qIJ"},
             const std::array<std::string,2>& noiseName = {"JrJC", "qJC"})
       : mtUnaryUpdate(name, errorName, stateName, noiseName),
         BrBC_(0,0,0),
         qBC_(1,0,0,0){
  }

  virtual ~PoseUpdate() {
  }

  void Eval(Vec3& JrJC_inn, Quat& qJC_inn,
            const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
            const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    JrJC_inn = ComputeExternalPosition(IrIB_cur, qIB_cur, IrIJ_cur, qIJ_cur)
        - meas_->JrJC_ + JrJC_noi;
    Quat dQ = dQ.exponentialMap(qJC_noi);
    qJC_inn = dQ * ComputeExternalAttitude(qIB_cur,qIJ_cur) * meas_->qJC_.inverted();
  }

  void JacCur(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    GetJacBlockCur<POS, POS>(J) = RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, ATT>(J) = -gSM(RotMat(qIJ_cur.inverted() * qIB_cur).rotate(BrBC_))
                                  *RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, IJP>(J) = -RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, IJA>(J) = RotMat(qIJ_cur).matrix().transpose()*gSM(Vec3(IrIB_cur
                                  - IrIJ_cur + qIB_cur.rotate(BrBC_)));
    GetJacBlockCur<ATT, ATT>(J) = RotMat(qIJ_cur.inverted()).matrix();
    GetJacBlockCur<ATT, IJA>(J) = -RotMat(qIJ_cur.inverted()).matrix();
  }

  void JacNoi(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    GetJacBlockNoi<POS, POS>(J) = Mat3::Identity();
    GetJacBlockNoi<ATT, ATT>(J) = Mat3::Identity();
  }

  void SetExtrinsics(Vec3 BrBC, Quat qBC){
    BrBC_  = BrBC;
    qBC_ = qBC;
  }

  void GetExtrinsics(Vec3* BrBC, Quat* qBC) const{
    *BrBC  =  BrBC_;
    *qBC = qBC_;
  }

  Vec3 ComputeExternalPosition(const Vec3& IrIB, const Quat& qIB,
                               const Vec3& IrIJ, const Quat& qIJ) const{
    return qIJ.inverseRotate(Vec3(IrIB - IrIJ + qIB.rotate(BrBC_)));
  }

  Vec3 ComputeExternalPosition(const ElementVectorBase& state) const{
    return ComputeExternalPosition(state.template GetValue<Vec3>(this->CurDefinition()->GetName(0)),
                                   state.template GetValue<Quat>(this->CurDefinition()->GetName(1)),
                                   state.template GetValue<Vec3>(this->CurDefinition()->GetName(2)),
                                   state.template GetValue<Quat>(this->CurDefinition()->GetName(3)));
  }

  Quat ComputeExternalAttitude(const Quat& qIB, const Quat& qIJ) const{
    return qIJ.inverted() * qIB * qBC_;
  };

  Quat ComputeExternalAttitude(const ElementVectorBase& state) const{
    return ComputeExternalAttitude(state.template GetValue<Quat>(this->CurDefinition()->GetName(1)),
                                   state.template GetValue<Quat>(this->CurDefinition()->GetName(3)));
  }

 protected:
  enum Elements {
    POS,
    ATT,
    IJP,
    IJA
  };
  Vec3 BrBC_;
  Quat qBC_;
};

}

#endif /* GIF_POSEUPDATE_HPP_ */
