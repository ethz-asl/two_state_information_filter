#ifndef GIF_RORMEAS_HPP_
#define GIF_RORMEAS_HPP_

#include "generalized_information_filter/common.h"

namespace GIF {

/*! \brief Rotational rate Measurement of Imu
 *         ElementVector that can be used to hold Rotational rate measurement.
 */
class RorMeas : public ElementVector {
 public:
  RorMeas(const Vec3& MwM = Vec3(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Vec3>({"MwM"}))),
        MwM_(ElementVector::GetValue<Vec3>("MwM")) {
    MwM_ = MwM;
  }
  Vec3& MwM_;
};

}
#endif /* GIF_RORMEAS_HPP_ */
