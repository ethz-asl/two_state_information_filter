#include "../include/generalized_information_filter/element-vector.h"

namespace GIF {

ElementVectorBase::ElementVectorBase(const std::shared_ptr<const ElementVectorDefinition>& def)
    : def_(def) {
}

ElementVectorBase::~ElementVectorBase() {}

bool ElementVectorBase::matchesDef(
    const std::shared_ptr<const ElementVectorDefinition>& def) const {
  if (getNumElement() != def->GetNumElements()) {
    return false;
  }
  for (int i = 0; i < getNumElement(); i++) {
    if (!def->GetElementDefinition(i)->isOfDef(getElement(i))) {
      return false;
    }
  }
  return true;
}

ElementVectorBase& ElementVectorBase::operator=(const ElementVectorBase& other) {
  for (int i = 0; i < getNumElement(); i++) {
    *getElement(i) = *other.getElement(i);
  }
  return *this;
}

void ElementVectorBase::print() const {
  for (int i = 0; i < getNumElement(); i++) {
    std::cout << getDef()->GetName(i) << ": ";
    getElement(i)->print();
  }
}

void ElementVectorBase::setIdentity() {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->setIdentity();
  }
}

void ElementVectorBase::setRandom(int& s) {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->setRandom(s);
  }
}

void ElementVectorBase::boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec,
                        const std::shared_ptr<ElementVectorBase>& out) const {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->boxplus(
        vec.block(getStart(i), 0, getElement(i)->getDim(), 1),
        out->getElement(i));
  }
}

void ElementVectorBase::boxminus(const std::shared_ptr<const ElementVectorBase>& ref,
                         Eigen::Ref<Eigen::VectorXd> vec) const {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->boxminus(
        ref->getElement(i),
        vec.block(getStart(i), 0, getElement(i)->getDim(), 1));
  }
}

std::shared_ptr<const ElementVectorDefinition> ElementVectorBase::getDef() const {
  return def_;
}

ElementVector::ElementVector(const std::shared_ptr<const ElementVectorDefinition>& def)
    : ElementVectorBase(def) {
  for (int i = 0; i < def_->GetNumElements(); i++) {
    elements_.push_back(def_->GetElementDefinition(i)->newElement());
  }
}

ElementVector::~ElementVector() {
}
ElementVector& ElementVector::operator=(const ElementVectorBase& other) {
  dynamic_cast<ElementVectorBase&>(*this) = other;
  return *this;
}

int ElementVector::getNumElement() const {
  return elements_.size();
}

ElementVectorWrapper::ElementVectorWrapper(const std::shared_ptr<const ElementVectorDefinition>& def,
                           const std::shared_ptr<const ElementVectorDefinition>& in)
    : ElementVectorBase(def),
      in_(in) {
  computeMap();
}

ElementVectorWrapper::~ElementVectorWrapper() {
}
ElementVectorWrapper& ElementVectorWrapper::operator=(const ElementVectorBase& other) {
  dynamic_cast<ElementVectorBase&>(*this) = other;
  return *this;
}

int ElementVectorWrapper::getNumElement() const {
  return indexMap_.size();
}

void ElementVectorWrapper::computeMap() {
  indexMap_.resize(def_->GetNumElements());
  for (int i = 0; i < def_->GetNumElements(); i++) {
    indexMap_[i] = in_->FindName(def_->GetName(i));
    assert(indexMap_.at(i) != -1);
  }
}

void ElementVectorWrapper::setState(const std::shared_ptr<ElementVectorBase>& state) {
  state->matchesDef(in_);
  state_ = state;
  constState_ = state;
}

void ElementVectorWrapper::setState(
    const std::shared_ptr<const ElementVectorBase>& state) const {
  state->matchesDef(in_);
  constState_ = state;
}

void ElementVectorWrapper::wrapJacobian(Eigen::Ref<MXD> out,
                                const Eigen::Ref<const MXD>& in,
                                int rowOffset) const {
  const int rows = in.rows();
  for (int i = 0; i < getNumElement(); ++i) {
    const int cols = getElement(i)->getDim();
    out.block(rowOffset, constState_->getStart(indexMap_.at(i)), rows, cols) =
        in.block(0, getStart(i), rows, cols);
  }
}

}
