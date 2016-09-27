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
  PoseMeas(const Vec3& JrJC = Vec3(0, 0, 0), const Quat& qCJ = Quat())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Quat>({ "JrJC", "qCJ" }))),
        JrJC_(ElementVector::GetValue<Vec3>("JrJC")),
        qCJ_(ElementVector::GetValue<Quat>("qCJ")) {
    JrJC_ = JrJC;
    qCJ_ = qCJ;
  }
  Vec3& JrJC_;
  Quat& qCJ_;
};

/*! \brief Pose Update
 *         Builds a residual between current estimate pose and an external measured pose. The state
 *         is parametrized by the position (IrIB) and attitude (qBI). Additionally, the offset
 *         between the inertial frame is coestimate (IrIJ and qJI). The extrinsics calibration of
 *         the pose measurement are assumed to be known  (BrBC and qCB).
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
  PoseUpdate(const std::array<std::string,2>& errorName = {"JrJC", "qCJ"},
             const std::array<std::string,4>& stateName = {"IrIB", "qBI", "IrIJ", "qJI"},
             const std::array<std::string,2>& noiseName = {"JrJC", "qCJ"})
       : mtUnaryUpdate(errorName, stateName, noiseName),
         BrBC_(0,0,0),
         qCB_(1,0,0,0){
  }

  virtual ~PoseUpdate() {
  }

  void Eval(Vec3& JrJC_inn, Quat& qCJ_inn,
            const Vec3& IrIB_cur, const Quat& qBI_cur, const Vec3& IrIJ_cur, const Quat& qJI_cur,
            const Vec3& JrJC_noi, const Vec3& qCJ_noi) const {
    JrJC_inn = qJI_cur.rotate(Vec3(IrIB_cur - IrIJ_cur + qBI_cur.inverseRotate(BrBC_)))
        - meas_->JrJC_ + JrJC_noi;
    Quat dQ = dQ.exponentialMap(qCJ_noi);
    qCJ_inn = dQ * qCB_* qBI_cur * qJI_cur.inverted() * meas_->qCJ_.inverted();
  }

  void JacCur(MatX& J,
              const Vec3& IrIB_cur, const Quat& qBI_cur, const Vec3& IrIJ_cur, const Quat& qJI_cur,
              const Vec3& JrJC_noi, const Vec3& qCJ_noi) const {
    J.setZero();
    GetJacBlockCur<POS, POS>(J) = RotMat(qJI_cur).matrix();
    GetJacBlockCur<POS, ATT>(J) = RotMat(qJI_cur * qBI_cur.inverted()).matrix()*gSM(BrBC_);
    GetJacBlockCur<POS, IJP>(J) = -RotMat(qJI_cur).matrix();
    GetJacBlockCur<POS, IJA>(J) = -gSM(qJI_cur.rotate(Vec3(IrIB_cur
                                  - IrIJ_cur + qBI_cur.inverseRotate(BrBC_))));
    GetJacBlockCur<ATT, ATT>(J) = RotMat(qCB_).matrix();
    GetJacBlockCur<ATT, IJA>(J) = -RotMat(qCB_* qBI_cur * qJI_cur.inverted()).matrix();
  }

  void JacNoi(MatX& J,
              const Vec3& IrIB_cur, const Quat& qBI_cur, const Vec3& IrIJ_cur, const Quat& qJI_cur,
              const Vec3& JrJC_noi, const Vec3& qCJ_noi) const {
    J.setZero();
    GetJacBlockNoi<POS, POS>(J) = Mat3::Identity();
    GetJacBlockNoi<ATT, ATT>(J) = Mat3::Identity();
  }

  void SetExtrinsics(Vec3 BrBC, Quat qCB){
    BrBC_  = BrBC;
    qCB_ = qCB;
  }

 protected:
  enum Elements {
    POS,
    ATT,
    IJP,
    IJA
  };
  Vec3 BrBC_;
  Quat qCB_;
};

}

#endif /* GIF_POSEUPDATE_HPP_ */
