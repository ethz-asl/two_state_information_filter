#ifndef GIF_ELEMENTVECTORDEFINITION_HPP_
#define GIF_ELEMENTVECTORDEFINITION_HPP_

#include <unordered_map>

#include "element-description.h"
#include "generalized_information_filter/common.h"

namespace GIF {

class ElementVectorBase;

class ElementVectorDefinition {
 public:
  ElementVectorDefinition();
  ~ElementVectorDefinition();
  bool operator ==(const std::shared_ptr<const ElementVectorDefinition>& other) const;
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
  void ExtendWithStateDefinition(const std::shared_ptr<const ElementVectorDefinition>& state_definition,
              const std::string& sub_name = "");

 protected:
  std::vector<std::pair<std::shared_ptr<const ElementDescriptionBase>, int>> element_definitions_;
  std::unordered_map<std::string, int> names_map_;
  int d_;
};


template<typename ... Ts>
struct TH_pack_dim;

template<typename ... Ts>
class ElementVectorPack : public ElementVectorDefinition{
 public:
  static constexpr int n_ = sizeof...(Ts);
  static constexpr int d_ = TH_pack_dim<Ts...>::d_;
  typedef std::tuple<Ts...> mtTuple;

  ElementVectorPack(const std::array<std::string,n_>& names);

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
int ElementVectorDefinition::GetStateDimension() const {
  return d_;
}

int ElementVectorDefinition::GetNumElements() const {
  return names_map_.size();
}

int ElementVectorDefinition::GetStartIndex(int outer_index) const {
  return element_definitions_.at(outer_index).second;
}

int ElementVectorDefinition::GetOuterIndex(int i) const {
  assert(i >= 0 && i < GetStateDimension());
  int outer_index = GetNumElements()-1;
  while (element_definitions_.at(outer_index).second > i) {
    --outer_index;
  }
  return outer_index;
}

int ElementVectorDefinition::GetInnerIndex(int i) const {
  return i - GetStartIndex(GetOuterIndex(i));
}

template<typename T>
int ElementVectorDefinition::AddElementDefinition(const std::string& name) {
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
ElementVectorPack<Ts...>::ElementVectorPack(const std::array<std::string,n_>& names){
  addElementsToDefinition(names);
}

template<typename ... Ts>
template<int i>
constexpr int ElementVectorPack<Ts...>::_GetStateDimension() {
  return ElementTraits<typename std::tuple_element<i,mtTuple>::type>::d_;
}

template<typename ... Ts>
template<int i, typename std::enable_if<(i>0)>::type*>
constexpr int ElementVectorPack<Ts...>::_GetStartIndex() {
  return _GetStartIndex<i-1>() + _GetStateDimension<i-1>();
}
template<typename ... Ts>
template<int i, typename std::enable_if<(i==0)>::type*>
constexpr int ElementVectorPack<Ts...>::_GetStartIndex() {
  return 0;
}

template<typename ... Ts>
template<int i, typename std::enable_if<(i<sizeof...(Ts))>::type*>
void ElementVectorPack<Ts...>::addElementsToDefinition(const std::array<std::string,n_>& names) {
  typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
  AddElementDefinition<mtElementType>(names.at(i));
  addElementsToDefinition<i+1>(names);
}
template<typename ... Ts>
template<int i, typename std::enable_if<(i==sizeof...(Ts))>::type*>
void ElementVectorPack<Ts...>::addElementsToDefinition(const std::array<std::string,n_>& names) {}

}

#endif /* GIF_ELEMENTVECTORDEFINITION_HPP_ */
