/*
 * State.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_STATE_HPP_
#define GIF_STATE_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Element.hpp"

namespace GIF{

class StateDefinition;

class State{
 public:
  State(StateDefinition* def): def_(def){};
  ~State(){
    for(auto e : elements_){
      delete e;
    }
  };
  State& operator=(const State& other){
    for(int i=0;i<elements_.size();i++){
      *elements_[i] = *other.elements_[i];
    }
    return *this;
  }
  void addElement(ElementBase* e){
    elements_.push_back(e);
  }
  ElementBase* getElement(int i){
    return elements_[i];
  }
  const ElementBase* getElement(int i) const{
    return elements_[i];
  }

 protected:
  std::vector<ElementBase*> elements_;
  StateDefinition* def_;
};

}

#endif /* GIF_STATE_HPP_ */
