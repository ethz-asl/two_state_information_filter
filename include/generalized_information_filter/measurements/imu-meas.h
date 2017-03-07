#ifndef GIF_IMUMEAS_HPP_
#define GIF_IMUMEAS_HPP_

#include "generalized_information_filter/common.h"

namespace GIF {

/*! \brief Imu Measurement
 *         ElementVector that can be used to hold Imu measurements.
 */
class ImuMeas : public ElementVector {
 public:
  ImuMeas(const Vec3& MwM = Vec3(0, 0, 0), const Vec3& MfM = Vec3(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Vec3>({"MwM", "MfM"}))),
        MwM_(ElementVector::GetValue<Vec3>("MwM")),
        MfM_(ElementVector::GetValue<Vec3>("MfM")) {
    MwM_ = MwM;
    MfM_ = MfM;
  }
  Vec3& MwM_;
  Vec3& MfM_;
};

}
#endif /* GIF_IMUMEAS_HPP_ */
