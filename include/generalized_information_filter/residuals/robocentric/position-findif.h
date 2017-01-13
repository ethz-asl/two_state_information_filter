#ifndef GIF_POSITIONFINDIF_HPP_
#define GIF_POSITIONFINDIF_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/measurement.h"

namespace GIF {

/*! \brief Robocentric velocity residual using finite differences.
 */
class PositionFindif
      : public BinaryResidual<ElementPack<Vec3>, ElementPack<Vec3, Vec3, Quat>,
                              ElementPack<Vec3>, ElementPack<Vec3>, EmptyMeas> {
 public:
  using mtResidual = BinaryResidual<ElementPack<Vec3>, ElementPack<Vec3, Vec3, Quat>,
                                    ElementPack<Vec3>, ElementPack<Vec3>, EmptyMeas>;
  PositionFindif(const std::string& name,
            const std::array<std::string,1>& innName = {"IrIM"},
            const std::array<std::string,3>& preName = {"IrIM", "MvM", "qIM"},
            const std::array<std::string,1>& curName = {"IrIM"},
            const std::array<std::string,1>& noiName = {"IrIM"})
    : mtResidual(name,innName,preName,curName,noiName,false,true,true){
  }
  virtual ~PositionFindif() {
  }
  enum Elements {POS, VEL, ATT};
  void Eval(Vec3& IrIM_inn, const Vec3& IrIM_pre, const Vec3& MvM_pre,
            const Quat& qIM_pre, const Vec3& IrIM_cur, const Vec3& IrIM_noi) const {
    IrIM_inn = IrIM_pre + dt_ * (qIM_pre.rotate(MvM_pre) + IrIM_noi / sqrt(dt_)) - IrIM_cur;
  }
  void JacPre(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre,
              const Quat& qIM_pre, const Vec3& IrIM_cur, const Vec3& IrIM_noi) const {
    J.setZero();
    this->template GetJacBlockPre<POS, POS>(J) = Mat3::Identity();
    this->template GetJacBlockPre<POS, VEL>(J) = dt_ * RotMat(qIM_pre).matrix();
    this->template GetJacBlockPre<POS, ATT>(J) = -dt_ * gSM(qIM_pre.rotate(MvM_pre));
  }
  void JacCur(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre,
              const Quat& qIM_pre, const Vec3& IrIM_cur, const Vec3& IrIM_noi) const {
    J.setZero();
    this->template GetJacBlockCur<POS, POS>(J) = -Mat3::Identity();
  }
  void JacNoi(MatX& J, const Vec3& IrIM_pre, const Vec3& MvM_pre,
              const Quat& qIM_pre, const Vec3& IrIM_cur, const Vec3& IrIM_noi) const {
    J.setZero();
    this->template GetJacBlockNoi<POS, POS>(J) = sqrt(dt_) * Mat3::Identity();
  }
};

}
#endif /* GIF_POSITIONFINDIF_HPP_ */
