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

namespace GIF{

class ElementDefinitionBase{
 public:
  ElementDefinitionBase(){};
  virtual ~ElementDefinitionBase(){};
  virtual std::shared_ptr<ElementBase> newElement() = 0;
  virtual bool isSame(const std::shared_ptr<const ElementDefinitionBase>& in) const = 0;
};

template<typename ElementType>
class ElementDefinition: public ElementDefinitionBase{
 public:
  ElementDefinition(){};
  ~ElementDefinition(){};
  std::shared_ptr<ElementBase> newElement(){
    return std::shared_ptr<Element<ElementType>>(new Element<ElementType>());
  }
  bool isSame(const std::shared_ptr<const ElementDefinitionBase>& in) const{
    return std::dynamic_pointer_cast<const ElementDefinition<ElementType>>(in).get() != nullptr;
  }
};

}

#endif /* GIF_ELEMENTDEFINITION_HPP_ */
