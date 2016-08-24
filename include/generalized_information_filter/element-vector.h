#ifndef GIF_ELEMENTVECTOR_HPP_
#define GIF_ELEMENTVECTOR_HPP_

#include "element-vector-definition.h"
#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"

namespace GIF {

class ElementVectorBase {
 public:
  ElementVectorBase(const std::shared_ptr<const ElementVectorDefinition>& def);
  virtual ~ElementVectorBase();
  bool matchesDef(const std::shared_ptr<const ElementVectorDefinition>& def) const;
  ElementVectorBase& operator=(const ElementVectorBase& other);
  virtual std::shared_ptr<ElementBase> getElement(int i) = 0;
  virtual std::shared_ptr<const ElementBase> getElement(int i) const = 0;
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
               const std::shared_ptr<ElementVectorBase>& out) const;
  void boxminus(const std::shared_ptr<const ElementVectorBase>& ref,
                Eigen::Ref<Eigen::VectorXd> vec) const;
  std::shared_ptr<const ElementVectorDefinition> getDef() const;

 protected:
  const std::shared_ptr<const ElementVectorDefinition> def_;
};

class ElementVector : public ElementVectorBase {
 public:
  ElementVector(const std::shared_ptr<const ElementVectorDefinition>& def);
  virtual ~ElementVector();
  ElementVector& operator=(const ElementVectorBase& other);
  int getNumElement() const;
  inline std::shared_ptr<ElementBase> getElement(int i);
  inline std::shared_ptr<const ElementBase> getElement(int i) const;

 protected:
  std::vector<std::shared_ptr<ElementBase>> elements_;
};

class ElementVectorWrapper : public ElementVectorBase {
 public:
  ElementVectorWrapper(const std::shared_ptr<const ElementVectorDefinition>& def,
               const std::shared_ptr<const ElementVectorDefinition>& in);
  ~ElementVectorWrapper();
  ElementVectorWrapper& operator=(const ElementVectorBase& other);
  int getNumElement() const;
  inline std::shared_ptr<ElementBase> getElement(int i);
  inline std::shared_ptr<const ElementBase> getElement(int i) const;
  void computeMap();
  void setState(const std::shared_ptr<ElementVectorBase>& state);
  void setState(const std::shared_ptr<const ElementVectorBase>& state) const;
  void wrapJacobian(Eigen::Ref<MXD> out, const Eigen::Ref<const MXD>& in,
                    int rowOffset = 0) const;

 protected:
  std::shared_ptr<ElementVectorBase> state_;
  mutable std::shared_ptr<const ElementVectorBase> constState_;
  const std::shared_ptr<const ElementVectorDefinition> in_;
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

std::shared_ptr<ElementBase> ElementVector::getElement(int i) {
  return elements_.at(i);
}

std::shared_ptr<const ElementBase> ElementVector::getElement(int i) const {
  return elements_.at(i);
}

std::shared_ptr<ElementBase> ElementVectorWrapper::getElement(int i) {
  return state_->getElement(indexMap_[i]);
}

std::shared_ptr<const ElementBase> ElementVectorWrapper::getElement(int i) const {
  return constState_->getElement(indexMap_[i]);
}

}
#endif /* GIF_ELEMENTVECTOR_HPP_ */
