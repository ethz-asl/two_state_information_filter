#ifndef GIF_ELEMENTVECTOR_HPP_
#define GIF_ELEMENTVECTOR_HPP_

#include "element-vector-definition.h"
#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"

namespace GIF {

class ElementVectorBase {
 public:
  ElementVectorBase(const SP<const ElementVectorDefinition>& def);
  virtual ~ElementVectorBase();
  bool matchesDef(const SP<const ElementVectorDefinition>& def) const;
  ElementVectorBase& operator=(const ElementVectorBase& other);
  virtual SP<ElementBase> getElement(int i) = 0;
  virtual SP<const ElementBase> getElement(int i) const = 0;
  template<typename T>
  inline T& getValue(int i);
  template<typename T>
  inline T& getValue(int i) const;
  template<typename T>
  T& getValue(const std::string& name);
  template<typename T>
  T& getValue(const std::string& name) const;
  inline int getDim() const;
  virtual int getNumElement() const = 0;
  inline int getStart(int i) const;
  inline int getOuter(int i) const;
  inline int getInner(int i) const;
  void print() const;
  void setIdentity();
  void setRandom(int& s);
  void boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec,
               const SP<ElementVectorBase>& out) const;
  void boxminus(const SP<const ElementVectorBase>& ref,
                Eigen::Ref<Eigen::VectorXd> vec) const;
  SP<const ElementVectorDefinition> getDef() const;

 protected:
  const SP<const ElementVectorDefinition> def_;
};

class ElementVector : public ElementVectorBase {
 public:
  ElementVector(const SP<const ElementVectorDefinition>& def);
  virtual ~ElementVector();
  ElementVector& operator=(const ElementVectorBase& other);
  int getNumElement() const;
  inline SP<ElementBase> getElement(int i);
  inline SP<const ElementBase> getElement(int i) const;

 protected:
  std::vector<SP<ElementBase>> elements_;
};

class ElementVectorWrapper : public ElementVectorBase {
 public:
  ElementVectorWrapper(const SP<const ElementVectorDefinition>& def,
               const SP<const ElementVectorDefinition>& in);
  ~ElementVectorWrapper();
  ElementVectorWrapper& operator=(const ElementVectorBase& other);
  int getNumElement() const;
  inline SP<ElementBase> getElement(int i);
  inline SP<const ElementBase> getElement(int i) const;
  void computeMap();
  void setState(const SP<ElementVectorBase>& state);
  void setState(const SP<const ElementVectorBase>& state) const;
  void wrapJacobian(Eigen::Ref<MXD> out, const Eigen::Ref<const MXD>& in,
                    int rowOffset = 0) const;

 protected:
  SP<ElementVectorBase> state_;
  mutable SP<const ElementVectorBase> constState_;
  const SP<const ElementVectorDefinition> in_;
  std::vector<int> indexMap_;
};

// ==================== Implementation ==================== //
template<typename T>
T& ElementVectorBase::getValue(int i) {
  return std::dynamic_pointer_cast < Element < T >> (getElement(i))->get();
}

template<typename T>
T& ElementVectorBase::getValue(int i) const {
  return std::dynamic_pointer_cast<const Element<T>>(getElement(i))->get();
}

template<typename T>
T& ElementVectorBase::getValue(const std::string& name) {
  assert(matchesDef(def_));
  int i = def_->FindName(name);
  assert(i != -1);
  return getValue<T>(i);
}

template<typename T>
T& ElementVectorBase::getValue(const std::string& name) const {
  assert(matchesDef(def_));
  int i = def_->FindName(name);
  assert(i != -1);
  return getValue<T>(i);
}

int ElementVectorBase::getDim() const {
  assert(matchesDef(def_));
  return def_->GetStateDimension();
}

int ElementVectorBase::getStart(int i) const {
  assert(matchesDef(def_));
  return def_->GetStartIndex(i);
}

int ElementVectorBase::getOuter(int i) const {
  assert(matchesDef(def_));
  return def_->GetOuterIndex(i);
}

int ElementVectorBase::getInner(int i) const {
  assert(matchesDef(def_));
  return def_->GetInnerIndex(i);
}

SP<ElementBase> ElementVector::getElement(int i) {
  return elements_.at(i);
}

SP<const ElementBase> ElementVector::getElement(int i) const {
  return elements_.at(i);
}

SP<ElementBase> ElementVectorWrapper::getElement(int i) {
  return state_->getElement(indexMap_[i]);
}

SP<const ElementBase> ElementVectorWrapper::getElement(int i) const {
  return constState_->getElement(indexMap_[i]);
}

}
#endif /* GIF_ELEMENTVECTOR_HPP_ */
