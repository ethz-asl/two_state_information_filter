/*
 * Element.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include "generalized_information_filter/common.hpp"

namespace GIF{

template<typename ElementType>
class ElementDefinition;

class ElementBase{
 public:
  ElementBase(){};
  virtual ~ElementBase(){};
  virtual ElementBase& operator=(const ElementBase& other) = 0;
};

template<typename ElementType>
class Element: public ElementBase{
 public:
  ElementType x_;
  Element(ElementDefinition<ElementType>* def): def_(def){};
  virtual ~Element(){};
  Element<ElementType>& operator=(const Element<ElementType>& other){
    x_ = other.x_;
    return *this;
  }
  virtual ElementBase& operator=(const ElementBase& other){
    *this = dynamic_cast<const Element<ElementType>&>(other);
    return *this;
  }
 protected:
  ElementDefinition<ElementType>* def_;
};

}

#endif /* ELEMENT_HPP_ */
