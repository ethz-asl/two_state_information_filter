#include "generalized_information_filter/StateDefinition.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

StateDefinition::StateDefinition(){
  d_ = 0;
};

StateDefinition::~StateDefinition(){};

bool StateDefinition::isSameDef(const std::shared_ptr<const StateDefinition>& in) const{
  if(getDim() != in->getDim()){
    return false;
  }
  for(auto name : namesMap_){
    int indexIn = in->findName(name.first);
    if(indexIn == -1){
      return false;
    }
    if(!getElementDefinition(name.second)->isSameDef(in->getElementDefinition(indexIn))){
      return false;
    }
  }
  return true;
}

std::string StateDefinition::getName(int i) const{
  for(auto e : namesMap_){
    if(e.second == i){
      return e.first;
    }
  }
  assert("Index not found in name map" == 0);
  return "";
}

int StateDefinition::findName(const std::string& name) const{
  auto query = namesMap_.find(name);
  return (query != namesMap_.end()) ? query->second : -1;
}

std::shared_ptr<const ElementDefinitionBase> StateDefinition::getElementDefinition(int i) const{
  return elementDefinitions_.at(i).first;
};

int StateDefinition::addElementDefinition(const std::string& name, const std::shared_ptr<const ElementDefinitionBase>& elementDefinition){
  int foundIndex = findName(name);
  if (foundIndex != -1){
    if(!getElementDefinition(foundIndex)->isSameDef(elementDefinition)){
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

void StateDefinition::extend(const std::shared_ptr<const StateDefinition>& stateDefinition, const std::string& name){
  for(auto nameEntry : stateDefinition->namesMap_){
    addElementDefinition(name + nameEntry.first, stateDefinition->getElementDefinition(nameEntry.second));
  }
}

}

