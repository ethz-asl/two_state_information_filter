#ifndef GIF_ELEMENTVECTOR_HPP_
#define GIF_ELEMENTVECTOR_HPP_

#include "generalized_information_filter/element-vector-definition.h"
#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"

namespace GIF {

/*! \brief Element Vector Base
 *         Base class providing the interface to a general Element Vector. Data
 *         may be stored inside or outside. Keeps a link to the employed
 *         definition.
 */
class ElementVectorBase {
 public:
  ElementVectorBase(const SP<const ElementVectorDefinition>& def);
  virtual ~ElementVectorBase();
  bool MatchesDefinition(const SP<const ElementVectorDefinition>& def) const;
  ElementVectorBase& operator=(const ElementVectorBase& other);
  virtual SP<ElementBase> GetElement(int i) = 0;
  virtual SP<const ElementBase> GetElement(int i) const = 0;
  template<typename T>
  inline T& GetValue(int i);
  template<typename T>
  inline T& GetValue(int i) const;
  template<typename T>
  T& GetValue(const std::string& name);
  template<typename T>
  T& GetValue(const std::string& name) const;
  inline int GetDimension() const;
  virtual int GetNumElement() const = 0;
  inline int GetStart(int i) const;
  inline int GetOuter(int i) const;
  inline int GetInner(int i) const;
  void Print() const;
  void SetIdentity();
  void SetRandom(int& s);
  void BoxPlus(const VecCRef<>& vec, const SP<ElementVectorBase>& out) const;
  void BoxMinus(const SP<const ElementVectorBase>& ref, VecRef<> vec) const;
  SP<const ElementVectorDefinition> GetDefinition() const;

 protected:
  const SP<const ElementVectorDefinition> def_;
};

/*! \brief Element Vector
 *         Holds the data as vector of element bases.
 */
class ElementVector : public ElementVectorBase {
 public:
  ElementVector(const SP<const ElementVectorDefinition>& def);
  virtual ~ElementVector();
  ElementVector& operator=(const ElementVectorBase& other);
  int GetNumElement() const;
  inline SP<ElementBase> GetElement(int i);
  inline SP<const ElementBase> GetElement(int i) const;

 protected:
  std::vector<SP<ElementBase>> elements_;
};

/*! \brief Element Vector Wrapper
 *         Holds a reference to an other element vector. This is used for
 *         mapping larger element vector into smaller ones which match a
 *         desired interface.
 */
class ElementVectorWrapper : public ElementVectorBase {
 public:
  ElementVectorWrapper(const SP<const ElementVectorDefinition>& def,
                       const SP<const ElementVectorDefinition>& in);
  ~ElementVectorWrapper();
  ElementVectorWrapper& operator=(const ElementVectorBase& other);
  int GetNumElement() const;
  inline SP<ElementBase> GetElement(int i);
  inline SP<const ElementBase> GetElement(int i) const;
  void computeMap();
  void setState(const SP<ElementVectorBase>& state);
  void setState(const SP<const ElementVectorBase>& state) const;
  void wrapJacobian(MatRef<> out, const MatCRef<>& in, int rowOffset = 0) const;

 protected:
  SP<ElementVectorBase> state_;
  mutable SP<const ElementVectorBase> constState_;
  const SP<const ElementVectorDefinition> inDef_;
  std::vector<int> indexMap_;
};

// ==================== Implementation ==================== //
template<typename T>
T& ElementVectorBase::GetValue(int i) {
  return std::dynamic_pointer_cast<Element<T>> (GetElement(i))->get();
}

template<typename T>
T& ElementVectorBase::GetValue(int i) const {
  return std::dynamic_pointer_cast<const Element<T>>(GetElement(i))->get();
}

template<typename T>
T& ElementVectorBase::GetValue(const std::string& name) {
  assert(MatchesDefinition(def_));
  int i = def_->FindName(name);
  assert(i != -1);
  return GetValue<T>(i);
}

template<typename T>
T& ElementVectorBase::GetValue(const std::string& name) const {
  assert(MatchesDefinition(def_));
  int i = def_->FindName(name);
  assert(i != -1);
  return GetValue<T>(i);
}

int ElementVectorBase::GetDimension() const {
  assert(MatchesDefinition(def_));
  return def_->GetStateDimension();
}

int ElementVectorBase::GetStart(int i) const {
  assert(MatchesDefinition(def_));
  return def_->GetStartIndex(i);
}

int ElementVectorBase::GetOuter(int i) const {
  assert(MatchesDefinition(def_));
  return def_->GetOuterIndex(i);
}

int ElementVectorBase::GetInner(int i) const {
  assert(MatchesDefinition(def_));
  return def_->GetInnerIndex(i);
}

SP<ElementBase> ElementVector::GetElement(int i) {
  return elements_.at(i);
}

SP<const ElementBase> ElementVector::GetElement(int i) const {
  return elements_.at(i);
}

SP<ElementBase> ElementVectorWrapper::GetElement(int i) {
  return state_->GetElement(indexMap_[i]);
}

SP<const ElementBase> ElementVectorWrapper::GetElement(int i) const {
  return constState_->GetElement(indexMap_[i]);
}

}
#endif /* GIF_ELEMENTVECTOR_HPP_ */
