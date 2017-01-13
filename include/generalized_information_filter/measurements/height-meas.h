#ifndef GIF_HEIGHTMEAS_HPP_
#define GIF_HEIGHTMEAS_HPP_

#include "generalized_information_filter/common.h"

namespace GIF {

/*! \brief Pose Measurement
 *         ElementVector that can be used to hold a generic pose measurements (position + attitude).
 */
class HeightMeas : public ElementVector {
 public:
  HeightMeas(const double& z = 0)
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(
            new ElementPack<double>({"z"}))),
        z_(ElementVector::GetValue<double>("z")) {
    z_ = z;
  }
  double& z_;
};

}

#endif /* GIF_HEIGHTMEAS_HPP_ */
