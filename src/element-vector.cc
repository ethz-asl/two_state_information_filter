#include "generalized_information_filter/element-vector.h"

namespace GIF {

ElementVectorBase::ElementVectorBase(
      const SP<const ElementVectorDefinition>& def): def_(def) {
}

ElementVectorBase::~ElementVectorBase() {}

bool ElementVectorBase::MatchesDefinition(
    const SP<const ElementVectorDefinition>& def) const {
  if (GetNumElement() != def->GetNumElements()) {
    return false;
  }
  for (int i = 0; i < GetNumElement(); i++) {
    if (!def->GetElementDefinition(i)->MatchesDescription(GetElement(i))) {
      return false;
    }
  }
  return true;
}

ElementVectorBase& ElementVectorBase::operator=(const ElementVectorBase& other){
  for (int i = 0; i < GetNumElement(); i++) {
    *GetElement(i) = *other.GetElement(i);
  }
  return *this;
}

void ElementVectorBase::Print() const {
  for (int i = 0; i < GetNumElement(); i++) {
    std::cout << GetDefinition()->GetName(i) << ": ";
    GetElement(i)->print();
  }
}

void ElementVectorBase::SetIdentity() {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->setIdentity();
  }
}

void ElementVectorBase::SetRandom(int& s) {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->setRandom(s);
  }
}

void ElementVectorBase::BoxPlus(const VecRC<>& vec,
                                const SP<ElementVectorBase>& out) const {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->boxplus(
        vec.block(GetStart(i), 0, GetElement(i)->getDim(), 1),
        out->GetElement(i));
  }
}

void ElementVectorBase::BoxMinus(const SP<const ElementVectorBase>& ref,
                         VecR<> vec) const {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->boxminus(
        ref->GetElement(i),
        vec.block(GetStart(i), 0, GetElement(i)->getDim(), 1));
  }
}

SP<const ElementVectorDefinition> ElementVectorBase::GetDefinition() const {
  return def_;
}

ElementVector::ElementVector(const SP<const ElementVectorDefinition>& def)
    : ElementVectorBase(def) {
  for (int i = 0; i < def_->GetNumElements(); i++) {
    elements_.push_back(def_->GetElementDefinition(i)->MakeElement());
  }
}

ElementVector::~ElementVector() {
}

ElementVector& ElementVector::operator=(const ElementVectorBase& other) {
  dynamic_cast<ElementVectorBase&>(*this) = other;
  return *this;
}

int ElementVector::GetNumElement() const {
  return elements_.size();
}

ElementVectorWrapper::ElementVectorWrapper(
      const SP<const ElementVectorDefinition>& def,
      const SP<const ElementVectorDefinition>& in): ElementVectorBase(def),
                                                    inDef_(in) {
  computeMap();
}

ElementVectorWrapper::~ElementVectorWrapper() {
}

ElementVectorWrapper& ElementVectorWrapper::operator=(
      const ElementVectorBase& other) {
  dynamic_cast<ElementVectorBase&>(*this) = other;
  return *this;
}

int ElementVectorWrapper::GetNumElement() const {
  return indexMap_.size();
}

void ElementVectorWrapper::computeMap() {
  indexMap_.resize(def_->GetNumElements());
  for (int i = 0; i < def_->GetNumElements(); i++) {
    indexMap_[i] = inDef_->FindName(def_->GetName(i));
    assert(indexMap_.at(i) != -1);
  }
}

void ElementVectorWrapper::setState(const SP<ElementVectorBase>& state) {
  state->MatchesDefinition(inDef_);
  state_ = state;
  constState_ = state;
}

void ElementVectorWrapper::setState(
      const SP<const ElementVectorBase>& state) const {
  state->MatchesDefinition(inDef_);
  constState_ = state;
}

void ElementVectorWrapper::wrapJacobian(MatR<> out,
                                        const MatRC<>& in,
                                        int rowOffset) const {
  const int rows = in.rows();
  for (int i = 0; i < GetNumElement(); ++i) {
    const int cols = GetElement(i)->getDim();
    out.block(rowOffset, constState_->GetStart(indexMap_.at(i)), rows, cols) =
        in.block(0, GetStart(i), rows, cols);
  }
}

}
