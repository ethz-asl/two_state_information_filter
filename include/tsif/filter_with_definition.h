#ifndef TSIF_FILTER_WITH_DEFINITION_H_
#define TSIF_FILTER_WITH_DEFINITION_H_

namespace tsif{

template<typename Definition>
class FilterWithDefinition : public Definition::FilterBase {
 public:

  using FD = Definition;
  using StateEnum = typename FD::StateEnum;
  using ParamEnum = typename FD::ParamEnum;
  using ResidualEnum = typename FD::ResidualEnum;

  explicit FilterWithDefinition(){}
  virtual ~FilterWithDefinition(){}

};

/*
#include <tsif/filter.h>
#include <tsif/residual_1.h>
#include <tsif/residual_1.h>
// etc.

Class FilterDefinitionExample {
  //filter state enum
  enum StateEnum : unsigned int {
    STATE_1 = 0,
    [...],
    NUM_STATES
  };
  //parameter (i.e. unoptimized state) enum
  //parameters need to have negative indices
  enum ParamEnum : int {
    PARAM_N = -N,
    [...],
    PARAM_1
  };
  //residual enum
  //this needs to match the order of the residual pack given to the tsif::Filter as template argument
  enum ResidualEnum : unsigned int {
    RESIDUAL_1 = 0,
    [...]
  };
  //base type for the filter
  //this needs to match the order of the residual enum defined above
  using FilterBase = Filter<[...]>;
};
*/

} // namespace tsif

#endif  // TSIF_FILTER_WITH_DEFINITION_H_
