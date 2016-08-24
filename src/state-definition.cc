#include "../include/generalized_information_filter/element-vector.h"
#include "../include/generalized_information_filter/element-vector-definition.h"

namespace GIF {

ElementVectorDefinition::ElementVectorDefinition() {
  d_ = 0;
}

ElementVectorDefinition::~ElementVectorDefinition() {}

bool ElementVectorDefinition::operator ==(
    const std::shared_ptr<const ElementVectorDefinition>& other) const {
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

std::string ElementVectorDefinition::GetName(int outer_index) const {
  for (auto e : names_map_) {
    if (e.second == outer_index) {
      return e.first;
    }
  }
  assert("Index not found in name map" == 0);
  return "";
}

int ElementVectorDefinition::FindName(const std::string& name) const {
  auto query = names_map_.find(name);
  return (query != names_map_.end()) ? query->second : -1;
}

std::shared_ptr<const ElementDescriptionBase> ElementVectorDefinition::GetElementDefinition(
    int outer_index) const {
  return element_definitions_.at(outer_index).first;
}

int ElementVectorDefinition::AddElementDefinition(
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

void ElementVectorDefinition::ExtendWithStateDefinition(
    const std::shared_ptr<const ElementVectorDefinition>& state_definition,
    const std::string& sub_name) {
  for (auto entry : state_definition->names_map_) {
    AddElementDefinition(
        sub_name + entry.first,
        state_definition->GetElementDefinition(entry.second));
  }
}

}

