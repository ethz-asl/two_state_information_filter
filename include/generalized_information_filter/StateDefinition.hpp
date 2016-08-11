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

class StateDefinition{ // TODO: shared_ptr
 public:
  StateDefinition(){};
  ~StateDefinition(){
    for(auto d : elementDefinitions_){
      delete d.second;
    }
  };
  State* newState() const{
    State* newState = new State();
    for(auto d : elementDefinitions_){
      newState->addElement(d.second->newElement());
    }
    return newState;
  }
  template<typename T>
  ElementDefinition<T>* addElementDefinition(const std::string& name){
    auto entry = elementDefinitions_.find(name);
    if (entry != elementDefinitions_.end()){
      return dynamic_cast<ElementDefinition<T>*>(entry->second);
    } else {
      ElementDefinition<T>* newDefinition = new ElementDefinition<T>();
      elementDefinitions_.insert(std::pair<std::string, ElementDefinitionBase*>(name,newDefinition));
      return newDefinition;
    }
  }

 protected:
  std::unordered_map<std::string, ElementDefinitionBase*> elementDefinitions_;
};

}

#endif /* GIF_STATEDEFINITION_HPP_ */
