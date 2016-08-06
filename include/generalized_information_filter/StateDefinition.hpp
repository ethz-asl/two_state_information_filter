/*
 * StateDefinition.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_STATEDEFINITION_HPP_
#define GIF_STATEDEFINITION_HPP_

#include <unordered_map>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/ElementDefinition.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

class StateDefinition{
 public:
  StateDefinition(){
    d_ = 0;
  };
  ~StateDefinition(){
    for(auto d : elementDefinitions_){
      delete d;
    }
  };
  State* newState(){
    State* newState = new State(this);
    for(auto d : elementDefinitions_){
      newState->addElement(d->newElement());
    }
    return newState;
  }
  template<typename T>
  ElementDefinition<T>* addElementDefinition(const std::string& name){
    if(namesMap_.count(name) == 0){
      ElementDefinition<T>* newDefinition = new ElementDefinition<T>(elementDefinitions_.size(),d_);
      d_ += newDefinition->getDim();
      elementDefinitions_.push_back(newDefinition);
      namesMap_[name] = newDefinition;
    }
    return dynamic_cast<ElementDefinition<T>*>(namesMap_.at(name));
  }
  int getDim(){
    return d_;
  }
  void print(const State* s) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i]->print(s->getElement(i));
    }
  }
  void init(State* s) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i]->init(s->getElement(i));
    }
  }
  void boxplus(const State* ref, const VXD& vec, State* out) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i]->boxplus(ref->getElement(i),vec.block(elementDefinitions_[i]->getIndex(),0,elementDefinitions_[i]->getDim(),1),out->getElement(i));
    }
  }
  void boxminus(const State* in, const State* ref, VXD& vec) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i]->boxminus(in->getElement(i),ref->getElement(i),vec.block(elementDefinitions_[i]->getIndex(),0,elementDefinitions_[i]->getDim(),1));
    }
  }

 protected:
  std::vector<ElementDefinitionBase*> elementDefinitions_;
  std::unordered_map<std::string, ElementDefinitionBase*> namesMap_;
  int d_;
};

}

#endif /* GIF_STATEDEFINITION_HPP_ */
