/*
 * ElementDefinition.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_ELEMENTDEFINITION_HPP_
#define GIF_ELEMENTDEFINITION_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Element.hpp"

namespace GIF {

class ElementDefinitionBase {
 public:
  ElementDefinitionBase() {}

  virtual ~ElementDefinitionBase() {}

  virtual std::shared_ptr<ElementBase> newElement() const = 0;
  virtual bool isSameDef(
      const std::shared_ptr<const ElementDefinitionBase>& in) const = 0;
  virtual bool isOfDef(const std::shared_ptr<const ElementBase>& in) const = 0;
  virtual int getDim() const = 0;
};

template<typename T>
class ElementDefinition : public ElementDefinitionBase {
 public:
  ElementDefinition() {
  }
  ;
  ~ElementDefinition() {
  }
  ;
  std::shared_ptr<ElementBase> newElement() const {
    return std::shared_ptr < Element < T >> (new Element<T>(this));
  }
  bool isSameDef(const std::shared_ptr<const ElementDefinitionBase>& in) const {
    return std::dynamic_pointer_cast<const ElementDefinition<T>>(in).get()
        != nullptr;
  }
  bool isOfDef(const std::shared_ptr<const ElementBase>& in) const {
    return std::dynamic_pointer_cast<const Element<T>>(in).get() != nullptr;
  }
  inline int getDim() const {
    return ElementTraits<T>::d_;
  }
};

}

#endif /* GIF_ELEMENTDEFINITION_HPP_ */
