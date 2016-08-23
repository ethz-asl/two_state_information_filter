#ifndef GIF_STATEDEFINITION_HPP_
#define GIF_STATEDEFINITION_HPP_

#include <unordered_map>

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element-definition.h"

namespace GIF {

class StateBase;

class StateDefinition {
 public:
  StateDefinition();
  ~StateDefinition();
  bool operator ==(const std::shared_ptr<const StateDefinition>& other) const;
  inline int GetStateDimension() const;
  inline int GetNumElements() const;
  inline int GetStartIndex(int outer_index) const;
  inline int GetOuterIndex(int i) const;
  inline int GetInnerIndex(int i) const;
  std::string GetName(int outer_index) const;
  int FindName(const std::string& name) const;
  std::shared_ptr<const ElementDefinitionBase> GetElementDefinition(
      int outer_index) const;
  int AddElementDefinition(
      const std::string& name,
      const std::shared_ptr<const ElementDefinitionBase>& new_element_definition);
  template<typename T>
  int AddElementDefinition(const std::string& name);
  void ExtendWithStateDefinition(const std::shared_ptr<const StateDefinition>& state_definition,
              const std::string& sub_name = "");

 protected:
  std::vector<std::pair<std::shared_ptr<const ElementDefinitionBase>, int>> element_definitions_;
  std::unordered_map<std::string, int> names_map_;
  int d_;
};

// ==================== Implementation ==================== //
int StateDefinition::GetStateDimension() const {
  return d_;
}

int StateDefinition::GetNumElements() const {
  return names_map_.size();
}

int StateDefinition::GetStartIndex(int outer_index) const {
  return element_definitions_.at(outer_index).second;
}

int StateDefinition::GetOuterIndex(int i) const {
  assert(i >= 0 && i < GetStateDimension());
  int outer_index = GetNumElements()-1;
  while (element_definitions_.at(outer_index).second > i) {
    --outer_index;
  }
  return outer_index;
}

int StateDefinition::GetInnerIndex(int i) const {
  return i - GetStartIndex(GetOuterIndex(i));
}

template<typename T>
int StateDefinition::AddElementDefinition(const std::string& name) {
  const std::shared_ptr<const ElementDefinitionBase> new_element_definition(
      new ElementDefinition<T>());
  AddElementDefinition(name, new_element_definition);
}

}

#endif /* GIF_STATEDEFINITION_HPP_ */
