#ifndef TSIF_FILTER_WITH_DEFINITION_H_
#define TSIF_FILTER_WITH_DEFINITION_H_

#include <tsif/filter.h>

namespace tsif{

//check if some type derives from a tsif
template <typename Derived>
struct is_tsif
{
  using D = typename std::remove_cv<Derived>::type;
  template <typename... Residuals>
  static std::true_type test(Filter<Residuals...>*);
  static std::false_type test(void*);
  using type = decltype(test(std::declval<D*>()));
};
template <typename Derived>
using is_tsif_t = typename is_tsif<Derived>::type;

template <typename ConcreteDefinition>
struct FilterDefinition {

  FilterDefinition() = delete;

  static_assert(!std::is_void<typename ConcreteDefinition::ParamEnum>::value,
                "[FilterDefinition] Filter definition must define ParamEnum.");
  using ParamEnum = typename ConcreteDefinition::ParamEnum;

  static_assert(!std::is_void<typename ConcreteDefinition::StateEnum>::value,
                "[FilterDefinition] Filter definition must define StateEnum.");
  using StateEnum = typename ConcreteDefinition::StateEnum;

  static_assert(!std::is_void<typename ConcreteDefinition::ResidualEnum>::value,
                "[FilterDefinition] Filter definition must define ResidualEnum.");
  using ResidualEnum = typename ConcreteDefinition::ResidualEnum;

  static_assert(tsif::is_tsif_t<typename ConcreteDefinition::FilterBase>::value,
                "[FilterDefinition] Filter definition must define FilterBase derived from a tsif::Filter");
  using FilterBase = typename ConcreteDefinition::FilterBase;

};

template<typename ConcreteDefinition>
class FilterWithDefinition : public ConcreteDefinition::FilterBase {
 public:

  using FD = FilterDefinition<ConcreteDefinition>;

  FilterWithDefinition() = default;
  ~FilterWithDefinition() override = default;

};

} // namespace tsif

#endif  // TSIF_FILTER_WITH_DEFINITION_H_
