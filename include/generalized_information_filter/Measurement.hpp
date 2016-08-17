/*
 * Measurement.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MEASUREMENT_HPP_
#define GIF_MEASUREMENT_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

class MeasurementBase: public State{
 public:
  MeasurementBase(const std::shared_ptr<const StateDefinition>& def): State(def){
    t0_ = 0;
    t1_ = 1;
  }
  virtual ~MeasurementBase(){};
  double t0_;
  double t1_;
};

}

#endif /* GIF_MEASUREMENT_HPP_ */
