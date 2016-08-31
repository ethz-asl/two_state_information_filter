#ifndef GIF_ELEMENTVECTOR_HPP_
#define GIF_ELEMENTVECTOR_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"
#include "generalized_information_filter/element-vector-definition.h"

namespace GIF {

/*! \brief Element Vector Base
 *         Base class providing the interface to a general Element Vector. Data
 *         may be stored inside or outside. Keeps a link to the employed
 *         definition.
 */
class ElementVectorBase {
 public:
  typedef std::shared_ptr<ElementVectorBase> Ptr;
  typedef std::shared_ptr<const ElementVectorBase> CPtr;
  ElementVectorBase(const ElementVectorDefinition::CPtr& def);
  virtual ~ElementVectorBase();
  bool MatchesDefinition(const ElementVectorDefinition::CPtr& def) const;
  ElementVectorBase& operator=(const ElementVectorBase& other);
  virtual ElementBase::Ptr GetElement(int i) = 0;
  virtual ElementBase::CPtr GetElement(int i) const = 0;
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
  void BoxPlus(const VecCRef<>& vec, const ElementVectorBase::Ptr& out) const;
  void BoxMinus(const ElementVectorBase::CPtr& ref, VecRef<> vec) const;
  ElementVectorDefinition::CPtr GetDefinition() const;

 protected:
  const ElementVectorDefinition::CPtr def_;
};

/*! \brief Element Vector
 *         Holds the data as vector of element bases.
 */
class ElementVector : public ElementVectorBase {
 public:
  ElementVector(const ElementVectorDefinition::CPtr& def);
  virtual ~ElementVector();
  ElementVector& operator=(const ElementVectorBase& other);
  int GetNumElement() const;
  inline ElementBase::Ptr GetElement(int i);
  inline ElementBase::CPtr GetElement(int i) const;

 protected:
  std::vector<ElementBase::Ptr> elements_;
};

/*! \brief Element Vector Wrapper
 *         Holds a reference to an other element vector. This is used for
 *         mapping larger element vector into smaller ones which match a
 *         desired interface.
 */
class ElementVectorWrapper : public ElementVectorBase {
 public:
  ElementVectorWrapper(const ElementVectorDefinition::CPtr& def,
                       const ElementVectorDefinition::CPtr& in);
  ~ElementVectorWrapper();
  ElementVectorWrapper& operator=(const ElementVectorBase& other);
  int GetNumElement() const;
  inline ElementBase::Ptr GetElement(int i);
  inline ElementBase::CPtr GetElement(int i) const;
  void ComputeMap();
  void SetVectorElement(const ElementVectorBase::Ptr& state);
  void SetVectorElement(const ElementVectorBase::CPtr& state) const;
  void EmbedJacobian(MatRef<> out, const MatCRef<>& in, int rowOffset = 0) const;

 protected:
  ElementVectorBase::Ptr state_;
  mutable ElementVectorBase::CPtr constState_;
  const ElementVectorDefinition::CPtr inDef_;
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

ElementBase::Ptr ElementVector::GetElement(int i) {
  return elements_.at(i);
}

ElementBase::CPtr ElementVector::GetElement(int i) const {
  return elements_.at(i);
}

ElementBase::Ptr ElementVectorWrapper::GetElement(int i) {
  return state_->GetElement(indexMap_[i]);
}

ElementBase::CPtr ElementVectorWrapper::GetElement(int i) const {
  return constState_->GetElement(indexMap_[i]);
}

}
#endif /* GIF_ELEMENTVECTOR_HPP_ */
