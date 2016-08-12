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
};

template<typename ElementType>
class ElementDefinition: public ElementDefinitionBase{
 public:
  ElementDefinition(){};
  ~ElementDefinition(){};
  std::shared_ptr<ElementBase> newElement(){
    return std::shared_ptr<Element<ElementType>>(new Element<ElementType>());
  }
};

}

#endif /* GIF_ELEMENTDEFINITION_HPP_ */
