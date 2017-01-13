#ifndef GIF_POSEMEAS_HPP_
#define GIF_POSEMEAS_HPP_

#include "generalized_information_filter/common.h"

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

}

#endif /* GIF_POSEMEAS_HPP_ */
