#ifndef GIF_POSEIBUPDATE_HPP_
#define GIF_POSEIBUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"
#include "generalized_information_filter/measurements/pose-meas.h"

namespace GIF {

/*! \brief Pose Update with time alignment
 *         Builds a residual between current estimate pose and an external measured pose. The state
 *         is parametrized by the position (IrIB) and attitude (qIB). Additionally, the offset
 *         between the inertial frame is coestimated (IrIJ and qIJ). The extrinsics calibration of
 *         the pose measurement is also coestimated  (BrBC and qBC).
 *
 *         Coordinate frames:
 *           B: Body
 *           C: Pose measurement sensor
 *           J: Pose measurement reference inertial frame
 *           I: Estimation reference inertial frame
 */
class PoseIBUpdate : public UnaryUpdate<ElementPack<Vec3, Quat>,
    ElementPack<Vec3, Quat, Vec3, Quat, Vec3, Quat>, ElementPack<Vec3, Vec3>, PoseMeas> {
 public:
  PoseIBUpdate(const std::string& name,
             const std::array<std::string,2>& errorName = {"JrJC", "qJC"},
             const std::array<std::string,6>& stateName = {"IrIB", "qIB", "IrIJ", "qIJ", "BrBC", "qBC"},
             const std::array<std::string,2>& noiseName = {"JrJC", "qJC"})
       : mtUnaryUpdate(name, errorName, stateName, noiseName){
  }

  virtual ~PoseIBUpdate() {
  }

  void Eval(Vec3& JrJC_inn, Quat& qJC_inn,
            const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
            const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    JrJC_inn = qIJ_cur.inverseRotate(Vec3(IrIB_cur - IrIJ_cur + qIB_cur.rotate(BrBC_cur)))
        - meas_->JrJC_ + JrJC_noi;
    Quat dQ = dQ.exponentialMap(qJC_noi);
    qJC_inn = dQ * qIJ_cur.inverted() * qIB_cur * qBC_cur * meas_->qJC_.inverted();
  }

  void JacCur(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    GetJacBlockCur<POS, POS>(J) = RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, ATT>(J) = -gSM(RotMat(qIJ_cur.inverted() * qIB_cur).rotate(BrBC_cur))
                                  *RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, IJP>(J) = -RotMat(qIJ_cur).matrix().transpose();
    GetJacBlockCur<POS, IJA>(J) = RotMat(qIJ_cur).matrix().transpose()*gSM(Vec3(IrIB_cur
                                  - IrIJ_cur + qIB_cur.rotate(BrBC_cur)));
    GetJacBlockCur<POS, BCP>(J) = RotMat(qIJ_cur.inverted()*qIB_cur).matrix();
    GetJacBlockCur<ATT, ATT>(J) = RotMat(qIJ_cur.inverted()).matrix();
    GetJacBlockCur<ATT, IJA>(J) = -RotMat(qIJ_cur.inverted()).matrix();
    GetJacBlockCur<ATT, BCA>(J) = RotMat(qIJ_cur.inverted() * qIB_cur).matrix();
  }

  void JacNoi(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    GetJacBlockNoi<POS, POS>(J) = Mat3::Identity();
    GetJacBlockNoi<ATT, ATT>(J) = Mat3::Identity();
  }

 protected:
  enum Elements {
    POS,
    ATT,
    IJP,
    IJA,
    BCP,
    BCA
  };
};

}

#endif /* GIF_POSEIBUPDATE_HPP_ */
