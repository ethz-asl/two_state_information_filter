#ifndef GIF_STATEDEFINITION_HPP_
#define GIF_STATEDEFINITION_HPP_

#include <unordered_map>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/ElementDefinition.hpp"

namespace GIF {

class StateBase;

class StateDefinition {
 public:
  StateDefinition();
  ~StateDefinition();
  bool isSameDef(const std::shared_ptr<const StateDefinition>& in) const;
  inline int getDim() const;
  inline int getNumElement() const;
  inline int getStart(int i) const;
  inline int getOuter(int i) const;  // TODO: make more efficient
  inline int getInner(int i) const;
  std::string getName(int i) const;
  int findName(const std::string& name) const;
  std::shared_ptr<const ElementDefinitionBase> getElementDefinition(
      int i) const;
  int addElementDefinition(
      const std::string& name,
      const std::shared_ptr<const ElementDefinitionBase>& elementDefinition);
  template<typename T>
  int addElementDefinition(const std::string& name);
  void extend(const std::shared_ptr<const StateDefinition>& stateDefinition,
              const std::string& name = "");

 protected:
  std::vector<std::pair<std::shared_ptr<const ElementDefinitionBase>, int>> elementDefinitions_;
  std::unordered_map<std::string, int> namesMap_;
  int d_;
};

// ==================== Implementation ==================== //
int StateDefinition::getDim() const {
  return d_;
}

int StateDefinition::getNumElement() const {
  return namesMap_.size();
}

int StateDefinition::getStart(int i) const {
  return elementDefinitions_.at(i).second;
}

int StateDefinition::getOuter(int i) const {  // TODO: make more efficient
  int j = 0;
  while (i >= getElementDefinition(j)->getDim()) {
    i -= getElementDefinition(j)->getDim();
    ++j;
  }
  return j;
}

int StateDefinition::getInner(int i) const {
  return i - getStart(getOuter(i));
}

template<typename T>
int StateDefinition::addElementDefinition(const std::string& name) {
  const std::shared_ptr<const ElementDefinitionBase> elementDefinition(
      new ElementDefinition<T>());
  addElementDefinition(name, elementDefinition);
}

}

#endif /* GIF_STATEDEFINITION_HPP_ */
