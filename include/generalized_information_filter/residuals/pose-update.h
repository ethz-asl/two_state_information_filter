#ifndef GIF_POSEUPDATE_HPP_
#define GIF_POSEUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"

namespace GIF {

class PoseMeas : public ElementVector {
 public:
  PoseMeas(const Vec3& pos = Vec3(0, 0, 0), const Quat& att = Quat())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Quat>({ "pos", "att" }))),
        pos_(ElementVector::GetValue<Vec3>("pos")),
        att_(ElementVector::GetValue<Quat>("att")) {
    pos_ = pos;
    att_ = att;
  }
  Vec3& pos_;
  Quat& att_;
};

class PoseUpdate : public UnaryUpdate<ElementPack<Vec3, Quat>,
    ElementPack<Vec3, Quat>, ElementPack<Vec3, Vec3>, PoseMeas> {
 public:
  PoseUpdate()
      : mtUnaryUpdate( { "pos", "att" }, { "pos", "att" }, { "pos", "att" }) {
    dt_ = 0.1;
  }

  virtual ~PoseUpdate() {
  }

  void eval(Vec3& posInn, Quat& attInn, const Vec3& posSta,
            const Quat& attSta, const Vec3& posNoi,
            const Vec3& attNoi) const {
    posInn = posSta - meas_->pos_ + posNoi;
    Quat dQ = dQ.exponentialMap(attNoi);
    attInn = dQ * attSta * meas_->att_.inverted();
  }
  void jacCur(MatX& J, const Vec3& posSta, const Quat& attSta,
                       const Vec3& posNoi, const Vec3& attNoi) const {
    J.setZero();
    setJacBlockCur<POS, POS>(J, Mat3::Identity());
    setJacBlockCur<ATT, ATT>(J, Mat3::Identity());
  }
  void jacNoi(MatX& J, const Vec3& posSta, const Quat& attSta,
                       const Vec3& posNoi, const Vec3& attNoi) const {
    J.setZero();
    setJacBlockNoi<POS, POS>(J, Mat3::Identity());
    setJacBlockNoi<ATT, ATT>(J, Mat3::Identity());
  }

 protected:
  double dt_;
  enum Elements {
    POS,
    ATT
  };
};

}

#endif /* GIF_POSEUPDATE_HPP_ */
