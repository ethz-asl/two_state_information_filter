#ifndef GIF_STATEDEFINITION_HPP_
#define GIF_STATEDEFINITION_HPP_

#include <unordered_map>

#include "element-description.h"
#include "generalized_information_filter/common.h"

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
  std::shared_ptr<const ElementDescriptionBase> GetElementDefinition(
      int outer_index) const;
  int AddElementDefinition(
      const std::string& name,
      const std::shared_ptr<const ElementDescriptionBase>& new_element_definition);
  template<typename T>
  int AddElementDefinition(const std::string& name);
  void ExtendWithStateDefinition(const std::shared_ptr<const StateDefinition>& state_definition,
              const std::string& sub_name = "");

 protected:
  std::vector<std::pair<std::shared_ptr<const ElementDescriptionBase>, int>> element_definitions_;
  std::unordered_map<std::string, int> names_map_;
  int d_;
};


template<typename ... Ts>
struct TH_pack_dim;

template<typename ... Ts>
class ElementPack : public StateDefinition{
 public:
  static constexpr int n_ = sizeof...(Ts);
  static constexpr int d_ = TH_pack_dim<Ts...>::d_;
  typedef std::tuple<Ts...> mtTuple;

  ElementPack(const std::array<std::string,n_>& names);

  template<int i>
  static constexpr int _GetStateDimension();

  template<int i, typename std::enable_if<(i>0)>::type* = nullptr>
  static constexpr int _GetStartIndex();
  template<int i, typename std::enable_if<(i==0)>::type* = nullptr>
  static constexpr int _GetStartIndex();

 protected:
  template<int i = 0, typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  void addElementsToDefinition(const std::array<std::string,n_>& names);
  template<int i = 0, typename std::enable_if<(i==sizeof...(Ts))>::type* = nullptr>
  void addElementsToDefinition(const std::array<std::string,n_>& names);

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
  const std::shared_ptr<const ElementDescriptionBase> new_element_definition(
      new ElementDescription<T>());
  AddElementDefinition(name, new_element_definition);
}

template<typename T, typename ... Ts>
struct TH_pack_dim<T, Ts...> {
  static constexpr int d_ = TH_pack_dim<Ts...>::d_ + ElementTraits<T>::d_;
};
template<>
struct TH_pack_dim<> {
  static constexpr int d_ = 0;
};

template<typename ... Ts>
ElementPack<Ts...>::ElementPack(const std::array<std::string,n_>& names){
  addElementsToDefinition(names);
}

template<typename ... Ts>
template<int i>
constexpr int ElementPack<Ts...>::_GetStateDimension() {
  return ElementTraits<typename std::tuple_element<i,mtTuple>::type>::d_;
}

template<typename ... Ts>
template<int i, typename std::enable_if<(i>0)>::type*>
constexpr int ElementPack<Ts...>::_GetStartIndex() {
  return _GetStartIndex<i-1>() + _GetStateDimension<i-1>();
}
template<typename ... Ts>
template<int i, typename std::enable_if<(i==0)>::type*>
constexpr int ElementPack<Ts...>::_GetStartIndex() {
  return 0;
}

template<typename ... Ts>
template<int i, typename std::enable_if<(i<sizeof...(Ts))>::type*>
void ElementPack<Ts...>::addElementsToDefinition(const std::array<std::string,n_>& names) {
  typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
  AddElementDefinition<mtElementType>(names.at(i));
  addElementsToDefinition<i+1>(names);
}
template<typename ... Ts>
template<int i, typename std::enable_if<(i==sizeof...(Ts))>::type*>
void ElementPack<Ts...>::addElementsToDefinition(const std::array<std::string,n_>& names) {}

}

#endif /* GIF_STATEDEFINITION_HPP_ */
