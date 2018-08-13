#ifndef TSIF_EXAMPLE_FILTER_H_
#define TSIF_EXAMPLE_FILTER_H_

// Filter definition
#include "tsif/example_filter_definition.h"

// Filter base
#include "tsif/filter_with_definition.h"

namespace tsif {

class ExampleFilter : public FilterWithDefinition<ExampleFilterDefinition> {
 public:

  using Base = FilterWithDefinition<ExampleFilterDefinition>;

  // Filter definition
  using Base::FD;

  ExampleFilter() = default;
  ~ExampleFilter() override = default;
  
};

} /* namespace tsif */

#endif  // TSIF_EXAMPLE_FILTER_H_