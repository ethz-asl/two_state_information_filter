#ifndef TSIF_EXAMPLE_FILTER_DEFINITION_H_
#define TSIF_EXAMPLE_FILTER_DEFINITION_H_

// Residuals
#include "tsif/residuals/random_walk.h"
// [...]

// Filter base
#include "tsif/filter.h"

namespace tsif {

struct ExampleFilterDefinition {
 public:
  
  // Definition contains only static members
  ExampleFilterDefinition() = delete;

  // Enum of uints to define the N+1 state variables
  enum StateEnum : unsigned int {
    STATE_0 = 0,
    STATE_1,
    // [...]
    STATE_N,
    // Number of states
    NUM_STATES
  };

  // Enum of negative ints to define the M+1 parameter variables
  enum ParamEnum : int {
    PARAM_M = -10, // PARAM_M = -(M+1), e.g. M+1 = 10
    PARAM_Mminus1,
    // [...]
    PARAM_0
  };

  // Enum of uints to define a filter made up of L+1 residuals
  enum ResidualEnum : unsigned int {
    RESIDUAL_0 = 0,
    RESIDUAL_1,
    // [...]
    RESIDUAL_L
  };

  // Base type for the filter. The order of template arguments must match the residual enum!
  using FilterBase = Filter<RandomWalk<Element<Vec3,STATE_0>>,
                            RandomWalk<Element<Vec3,STATE_1>>,
                            RandomWalk<Element<Vec3,STATE_N>>>;

};

} /* namespace tsif */

#endif  // TSIF_EXAMPLE_FILTER_DEFINITION_H_