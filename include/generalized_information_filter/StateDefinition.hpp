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

namespace GIF{

class StateDefinition{
 public:
  StateDefinition(){
    d_ = 0;
  };
  ~StateDefinition(){};
  bool isSame(const std::shared_ptr<const StateDefinition>& in) const{
    if(getDim() != in->getDim()){
      return false;
    }
    for(auto name : namesMap_){
      int indexIn = in->findName(name.first);
      if(indexIn == -1){
        return false;
      }
      if(!getElementDefinition(name.second)->isSame(in->getElementDefinition(indexIn))){
        return false;
      }
    }
    return true;
  }
  inline int getDim() const{
    return d_;
  }
  inline int getNumElement() const{
    return namesMap_.size();
  }
  inline int getStart(int i) const{
    return elementDefinitions_.at(i).second;
  }
  inline int getOuter(int i) const{ // TODO: make more efficient
    int j=0;
    while(i >= getElementDefinition(j)->getDim()){
      i -= getElementDefinition(j)->getDim();
      ++j;
    }
    return j;
  }
  inline int getInner(int i) const{
    return i-getStart(getOuter(i));
  }
  std::string getName(int i) const{
    for(auto e : namesMap_){
      if(e.second == i){
        return e.first;
      }
    }
    assert("Index not found in name map" == 0);
    return "";
  }
  int findName(const std::string& name) const{
    auto query = namesMap_.find(name);
    return (query != namesMap_.end()) ? query->second : -1;
  }
  const std::shared_ptr<const ElementDefinitionBase>& getElementDefinition(int i) const{
    return elementDefinitions_.at(i).first;
  };
  int addElementDefinition(const std::string& name, const std::shared_ptr<const ElementDefinitionBase>& elementDefinition){
    int foundIndex = findName(name);
    if (foundIndex != -1){
      if(!getElementDefinition(foundIndex)->isSame(elementDefinition)){
        assert("ERROR: invalid extention of state definition" == 0);
      }
      return foundIndex;
    } else {
      elementDefinitions_.push_back(std::pair<std::shared_ptr<const ElementDefinitionBase>,int>(elementDefinition,d_));
      d_ += elementDefinition->getDim();
      namesMap_.insert(std::pair<std::string, int>(name,elementDefinitions_.size()-1));
      return elementDefinitions_.size()-1;
    }
  }
  template<typename T>
  int addElementDefinition(const std::string& name){
    const std::shared_ptr<const ElementDefinitionBase> elementDefinition(new ElementDefinition<T>());
    addElementDefinition(name,elementDefinition);
  }
  void extend(const std::shared_ptr<const StateDefinition>& stateDefinition, const std::string& name = ""){
    for(auto nameEntry : stateDefinition->namesMap_){
      addElementDefinition(name + nameEntry.first, stateDefinition->getElementDefinition(nameEntry.second));
    }
  }

 protected:
  std::vector<std::pair<std::shared_ptr<const ElementDefinitionBase>,int>> elementDefinitions_;
  std::unordered_map<std::string, int> namesMap_;
  int d_;
};

}

#endif /* GIF_STATEDEFINITION_HPP_ */
