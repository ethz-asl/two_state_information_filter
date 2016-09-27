#ifndef GIF_ELEMENTVECTORDEFINITION_HPP_
#define GIF_ELEMENTVECTORDEFINITION_HPP_

#include <unordered_map>

#include "generalized_information_filter/element-description.h"
#include "generalized_information_filter/common.h"

namespace GIF {

class ElementVectorBase;

/*! \brief Element Vector Definition
 *         Defines the structure of an element vector. The descriptions of the
 *         elements are stored as vector element_definions_. mames_map_ stores
 *         the identifiers (strings) of the elements and provides the map to
 *         the outer index.
 */
class ElementVectorDefinition {
 public:
  typedef std::shared_ptr<ElementVectorDefinition> Ptr;
  typedef std::shared_ptr<const ElementVectorDefinition> CPtr;
  ElementVectorDefinition();
  virtual ~ElementVectorDefinition();
  bool MatchesDefinition(const ElementVectorDefinition& other) const;
  bool MatchesDefinition(const ElementVectorBase& other) const;
  inline int GetDim() const;
  inline int GetNumElements() const;
  inline int GetStart(int outer_index) const;
  inline int GetOuter(int i) const;
  inline int GetInner(int i) const;
  std::string GetName(int outer_index) const;
  int FindName(const std::string& name) const;
  const ElementDescriptionBase::CPtr& GetElementDescription(int outer_index) const;
  int AddElement(const std::string& name,
                 const ElementDescriptionBase::CPtr& new_element_definition);
  template<typename T>
  int AddElement(const std::string& name);
  void ExtendWithOtherElementVectorDefinition(const ElementVectorDefinition& other);

 protected:
  std::vector<ElementDescriptionBase::CPtr> descriptions_;
  std::vector<int> start_indices_;
  std::unordered_map<std::string, int> names_map_;
  int dim_;
};

/*! \brief Template Helper class for computing the dimension of an ElementPack
 */
template<typename ... Ts>
struct TH_pack_dim;

/*! \brief Element Pack
 *         This is a compile-time version of ElementVectorDefinition and
 *         and provides extra functionality
 */
template<typename ... Ts>
class ElementPack : public ElementVectorDefinition{
 public:
  static constexpr int n_ = sizeof...(Ts);
  static constexpr int kDim = TH_pack_dim<Ts...>::kDim;
  typedef std::tuple<Ts...> Tuple;

  ElementPack(const std::array<std::string,n_>& names);

  template<int i>
  static constexpr int _GetStateDimension();

  template<int i, typename std::enable_if<(i>0)>::type* = nullptr>
  static constexpr int _GetStartIndex();
  template<int i, typename std::enable_if<(i==0)>::type* = nullptr>
  static constexpr int _GetStartIndex();

 protected:
  template<int i = 0,
      typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  void addElementsToDefinition(const std::array<std::string,n_>& names);
  template<int i = 0,
        typename std::enable_if<(i==sizeof...(Ts))>::type* = nullptr>
  void addElementsToDefinition(const std::array<std::string,n_>& names);

};

// ==================== Implementation ==================== //
int ElementVectorDefinition::GetDim() const {
  return dim_;
}

int ElementVectorDefinition::GetNumElements() const {
  return descriptions_.size();
}

int ElementVectorDefinition::GetStart(int outer_index) const {
  return start_indices_.at(outer_index);
}

int ElementVectorDefinition::GetOuter(int i) const {
  DLOG_IF(FATAL, i < 0 | i >= GetDim()) << "Index out of range";
  int outer_index = GetNumElements()-1;
  while (GetStart(outer_index) > i) {
    --outer_index;
  }
  return outer_index;
}

int ElementVectorDefinition::GetInner(int i) const {
  return i - GetStart(GetOuter(i));
}

template<typename T>
int ElementVectorDefinition::AddElement(const std::string& name) {
  AddElement(name, std::make_shared<ElementDescription<T>>());
}

template<typename T, typename ... Ts>
struct TH_pack_dim<T, Ts...> {
  static constexpr int kDim = TH_pack_dim<Ts...>::kDim + ElementTraits<T>::kDim;
};
template<>
struct TH_pack_dim<> {
  static constexpr int kDim = 0;
};

template<typename ... Ts>
ElementPack<Ts...>::ElementPack(const std::array<std::string,n_>& names){
  addElementsToDefinition(names);
}

template<typename ... Ts>
template<int i>
constexpr int ElementPack<Ts...>::_GetStateDimension() {
  typedef typename std::tuple_element<i,Tuple>::type mtElementType;
  return ElementTraits<mtElementType>::kDim;
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
void ElementPack<Ts...>::addElementsToDefinition(
      const std::array<std::string,n_>& names) {
  typedef typename std::tuple_element<i,Tuple>::type mtElementType;
  AddElement<mtElementType>(names.at(i));
  addElementsToDefinition<i+1>(names);
}
template<typename ... Ts>
template<int i, typename std::enable_if<(i==sizeof...(Ts))>::type*>
void ElementPack<Ts...>::addElementsToDefinition(
      const std::array<std::string,n_>& names) {}

}

#endif /* GIF_ELEMENTVECTORDEFINITION_HPP_ */
