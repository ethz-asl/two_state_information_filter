#ifndef GIF_ACCMEAS_HPP_
#define GIF_ACCMEAS_HPP_

#include "generalized_information_filter/common.h"

namespace GIF {

/*! \brief Accelerometer Measurement of Imu
 *         ElementVector that can be used to hold acceleromter measurement.
 */
class AccMeas : public ElementVector {
 public:
  AccMeas(const Vec3& MfM = Vec3(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<Vec3, Vec3>({"MfM"}))),
        MfM_(ElementVector::GetValue<Vec3>("MfM")) {
    MfM_ = MfM;
  }
  Vec3& MfM_;
};

}
#endif /* GIF_ACCMEAS_HPP_ */
