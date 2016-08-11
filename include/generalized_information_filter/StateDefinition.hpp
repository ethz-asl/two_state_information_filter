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

class StateDefinition{ // TODO: hshared ptrs
 public:
  StateDefinition(){
    d_ = 0;
  };
  ~StateDefinition(){
    for(auto d : elementDefinitions_){
      delete d.first;
    }
  };
  State* newState() const{
    State* newState = new State(this);
    for(auto d : elementDefinitions_){
      newState->addElement(d.first->newElement());
    }
    return newState;
  }
  template<typename T>
  ElementDefinition<T>* addElementDefinition(const std::string& name){
    auto q = namesMap_.find(name);
    if (q != namesMap_.end()){
      return dynamic_cast<ElementDefinition<T>*>(elementDefinitions_.at(q->second).first);
    } else {
      namesMap_[name] = elementDefinitions_.size();
      ElementDefinition<T>* newDefinition = new ElementDefinition<T>();
      elementDefinitions_.push_back(std::pair<ElementDefinitionBase*,int>(newDefinition,d_));
      d_ += newDefinition->getDim();
      return newDefinition;
    }
  }
  int getDim() const{
    return d_;
  }
  void print(const State* s) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i].first->print(s->getElement(i));
    }
  }
  void init(State* s) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i].first->init(s->getElement(i));
    }
  }
  void boxplus(const State* ref, const VXD& vec, State* out) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i].first->boxplus(ref->getElement(i),vec.block(elementDefinitions_[i].second,0,elementDefinitions_[i].first->getDim(),1),out->getElement(i));
    }
  }
  void boxminus(const State* in, const State* ref, VXD& vec) const{
    for(int i=0;i<elementDefinitions_.size();i++){
      elementDefinitions_[i].first->boxminus(in->getElement(i),ref->getElement(i),vec.block(elementDefinitions_[i].second,0,elementDefinitions_[i].first->getDim(),1));
    }
  }

 protected:
  std::vector<std::pair<ElementDefinitionBase*,int>> elementDefinitions_;
  std::unordered_map<std::string, int> namesMap_;
  int d_;
};

}

#endif /* GIF_STATEDEFINITION_HPP_ */
