#include "generalized_information_filter/state-definition.h"
#include "generalized_information_filter/state.h"

namespace GIF {

StateDefinition::StateDefinition() {
  d_ = 0;
}

StateDefinition::~StateDefinition() {}

bool StateDefinition::operator ==(
    const std::shared_ptr<const StateDefinition>& other) const {
  if (GetStateDimension() != other->GetStateDimension()) {
    return false;
  }
  for (auto name : names_map_) {
    int other_outer_index = other->FindName(name.first);
    if (other_outer_index == -1) {
      return false;
    }
    if (!GetElementDefinition(name.second)->isSameDef(
        other->GetElementDefinition(other_outer_index))) {
      return false;
    }
  }
  return true;
}

std::string StateDefinition::GetName(int outer_index) const {
  for (auto e : names_map_) {
    if (e.second == outer_index) {
      return e.first;
    }
  }
  assert("Index not found in name map" == 0);
  return "";
}

int StateDefinition::FindName(const std::string& name) const {
  auto query = names_map_.find(name);
  return (query != names_map_.end()) ? query->second : -1;
}

std::shared_ptr<const ElementDescriptionBase> StateDefinition::GetElementDefinition(
    int outer_index) const {
  return element_definitions_.at(outer_index).first;
}

int StateDefinition::AddElementDefinition(
    const std::string& name,
    const std::shared_ptr<const ElementDescriptionBase>& new_element_definition) {
  int outer_index = FindName(name);
  if (outer_index != -1) {
    if (!GetElementDefinition(outer_index)->isSameDef(new_element_definition)) {
      assert("ERROR: invalid extention of state definition" == 0);
    }
    return outer_index;
  } else {
    element_definitions_.push_back(
        std::pair<std::shared_ptr<const ElementDescriptionBase>, int>(
            new_element_definition, d_));
    d_ += new_element_definition->getDim();
    names_map_.insert(
        std::pair<std::string, int>(name, element_definitions_.size() - 1));
    return element_definitions_.size() - 1;
  }
}

void StateDefinition::ExtendWithStateDefinition(
    const std::shared_ptr<const StateDefinition>& state_definition,
    const std::string& sub_name) {
  for (auto entry : state_definition->names_map_) {
    AddElementDefinition(
        sub_name + entry.first,
        state_definition->GetElementDefinition(entry.second));
  }
}

}

