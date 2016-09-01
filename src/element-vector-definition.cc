#include "generalized_information_filter/element-vector.h"
#include "generalized_information_filter/element-vector-definition.h"

namespace GIF {

ElementVectorDefinition::ElementVectorDefinition() {
  dim_ = 0;
}

ElementVectorDefinition::~ElementVectorDefinition() {}

bool ElementVectorDefinition::MatchesDefinition(const ElementVectorDefinition& other) const {
  if (GetDim() != other.GetDim()) {
    return false;
  }
  for (const auto& entry : names_map_) {
    int other_outer_index = other.FindName(entry.first);
    if (other_outer_index != entry.second) {
      return false;
    }
    if (!GetElementDescription(entry.second).MatchesDescription(
        other.GetElementDescription(other_outer_index))) {
      return false;
    }
  }
  return true;
}

bool ElementVectorDefinition::MatchesDefinition(
    const ElementVectorBase& other) const {
  return MatchesDefinition(*other.GetDefinition());
}

std::string ElementVectorDefinition::GetName(int outer_index) const { // Slow
  for (const auto& e : names_map_) {
    if (e.second == outer_index) {
      return e.first;
    }
  }
  assert("Index not found in name map" == 0);
  return "";
}

int ElementVectorDefinition::FindName(const std::string& name) const {
  const auto& query = names_map_.find(name);
  return (query != names_map_.end()) ? query->second : -1;
}

const ElementDescriptionBase& ElementVectorDefinition::GetElementDescription(int outer_index) const {
  return *descriptions_.at(outer_index);
}

int ElementVectorDefinition::AddElement(const std::string& name,
                                        const ElementDescriptionBase& description) {
  int outer_index = FindName(name);
  if (outer_index != -1) {
    if (!GetElementDescription(outer_index).MatchesDescription(description)) {
      assert("ERROR: invalid extension of state definition" == 0);
    }
    return outer_index;
  } else {
    descriptions_.push_back(description.Copy());
    start_indices_.push_back(dim_);
    dim_ += description.GetDimension();
    names_map_.insert(std::pair<std::string, int>(name, GetNumElements() - 1));
    return GetNumElements() - 1;
  }
}

void ElementVectorDefinition::ExtendWithOtherElementVectorDefinition(
      const ElementVectorDefinition& other) {
  for (const auto& entry : other.names_map_) {
    AddElement(entry.first,other.GetElementDescription(entry.second));
  }
}

}

