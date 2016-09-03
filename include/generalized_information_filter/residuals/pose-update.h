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
  PoseMeas(const Vec3& IrIB = Vec3(0, 0, 0), const Quat& qIB = Quat())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Quat>({ "IrIB", "qIB" }))),
        IrIB_(ElementVector::GetValue<Vec3>("IrIB")),
        qIB_(ElementVector::GetValue<Quat>("qIB")) {
    IrIB_ = IrIB;
    qIB_ = qIB;
  }
  Vec3& IrIB_;
  Quat& qIB_;
};

/*! \brief Pose Update
 *         Builds a residual between current state pose and measured pose.
 *         Coordinate frames: I (world), B (IMU).
 *         Chosen state parametrisation: Position (IrIB), Attitude (qIB).
 */
class PoseUpdate : public UnaryUpdate<ElementPack<Vec3, Quat>,
    ElementPack<Vec3, Quat>, ElementPack<Vec3, Vec3>, PoseMeas> {
 public:
  PoseUpdate()
      : mtUnaryUpdate( { "IrIB", "qIB" }, { "IrIB", "qIB" }, { "IrIB", "qIB" }) {
  }

  virtual ~PoseUpdate() {
  }

  void Eval(Vec3& pos_inn, Quat& att_inn, const Vec3& IrIB_cur,
            const Quat& qIB_cur, const Vec3& pos_noi,
            const Vec3& att_noi) const {
    pos_inn = IrIB_cur - meas_->IrIB_ + pos_noi;
    Quat dQ = dQ.exponentialMap(att_noi);
    att_inn = dQ * qIB_cur * meas_->qIB_.inverted();
  }
  void JacCur(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
                       const Vec3& pos_noi, const Vec3& att_noi) const {
    J.setZero();
    GetJacBlockCur<POS, POS>(J) = Mat3::Identity();
    GetJacBlockCur<ATT, ATT>(J) = Mat3::Identity();
  }
  void JacNoi(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
                       const Vec3& pos_noi, const Vec3& att_noi) const {
    J.setZero();
    GetJacBlockNoi<POS, POS>(J) = Mat3::Identity();
    GetJacBlockNoi<ATT, ATT>(J) = Mat3::Identity();
  }

 protected:
  enum Elements {
    POS,
    ATT
  };
};

}

#endif /* GIF_POSEUPDATE_HPP_ */
