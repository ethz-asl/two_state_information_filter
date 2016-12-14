#ifndef GIF_ROBCENROTRATEFINDIFRESIDUAL_HPP_
#define GIF_ROBCENROTRATEFINDIFRESIDUAL_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"

namespace GIF {

/*! \brief Robocentric rotational rate residual using finite differences.
 */
class RobcenRotrateFindifResidual
      : public BinaryResidual<ElementPack<Vec3>, ElementPack<Quat>,
                              ElementPack<Quat, Vec3>, ElementPack<Vec3>, EmptyMeas> {
 public:
  using mtResidual = BinaryResidual<ElementPack<Vec3>, ElementPack<Quat>,
                                    ElementPack<Quat, Vec3>, ElementPack<Vec3>, EmptyMeas>;
  RobcenRotrateFindifResidual(const std::string& name,
            const std::array<std::string,1>& innName = {"qIM"},
            const std::array<std::string,1>& preName = {"qIM"},
            const std::array<std::string,2>& curName = {"qIM", "MwM"},
            const std::array<std::string,1>& noiName = {"qIM"})
    : mtResidual(name,innName,preName,curName,noiName,false,true,true){
  }
  virtual ~RobcenRotrateFindifResidual() {
  }
  enum Elements {ATT, ROR};
  void Eval(Vec3& qIM_inn, const Quat& qIM_pre,
            const Quat& qIM_cur, const Vec3& MwM_cur, const Vec3& qIM_noi) const {
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur + qIM_noi * sqrt(dt_));
    qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
  }
  void JacPre(MatX& J, const Quat& qIM_pre,
              const Quat& qIM_cur, const Vec3& MwM_cur, const Vec3& qIM_noi) const {
    J.setZero();
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockPre<ATT, ATT>(J) = GammaMat(qIM_inn).inverse();
  }
  void JacCur(MatX& J, const Quat& qIM_pre,
              const Quat& qIM_cur, const Vec3& MwM_cur, const Vec3& qIM_noi) const {
    J.setZero();
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockCur<ATT, ATT>(J) = -GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre * dQ * qIM_cur.inverted()).matrix();
    this->template GetJacBlockCur<ATT, ROR>(J) = dt_ * GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cur);
  }
  void JacNoi(MatX& J, const Quat& qIM_pre,
              const Quat& qIM_cur, const Vec3& MwM_cur, const Vec3& qIM_noi) const {
    J.setZero();
    Quat dQ = dQ.exponentialMap(dt_ * MwM_cur);
    Vec3 qIM_inn = (qIM_pre * dQ).boxMinus(qIM_cur);
    this->template GetJacBlockNoi<ATT, ATT>(J) = sqrt(dt_) * GammaMat(qIM_inn).inverse() *
        RotMat(qIM_pre).matrix() * GammaMat(dt_ * MwM_cur);
  }
};

}
#endif /* GIF_ROBCENROTRATEFINDIFRESIDUAL_HPP_ */
