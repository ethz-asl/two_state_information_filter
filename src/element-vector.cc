#include "generalized_information_filter/element-vector.h"

namespace GIF {

ElementVectorBase::ElementVectorBase(
      const ElementVectorDefinition::CPtr& def): def_(def) {
}

ElementVectorBase::~ElementVectorBase() {}

bool ElementVectorBase::MatchesDefinition(
    const ElementVectorDefinition& def) const {
  if (GetNumElement() != def.GetNumElements()) {
    return false;
  }
  for (int i = 0; i < GetNumElement(); i++) {
    if (!def.GetElementDescription(i).MatchesDescription(*GetElement(i))) {
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
    GetElement(i)->Print();
  }
}

void ElementVectorBase::SetIdentity() {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->SetIdentity();
  }
}

void ElementVectorBase::SetRandom() {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->SetRandom();
  }
}

void ElementVectorBase::BoxPlus(const VecCRefX& vec, ElementVectorBase* out) const {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->Boxplus(
        vec.block(GetStart(i), 0, GetElement(i)->GetDim(), 1),
        out->GetElement(i));
  }
}

void ElementVectorBase::BoxMinus(const ElementVectorBase& ref,
                         VecRefX vec) const {
  for (int i = 0; i < GetNumElement(); i++) {
    GetElement(i)->Boxminus(
        *ref.GetElement(i),
        vec.block(GetStart(i), 0, GetElement(i)->GetDim(), 1));
  }
}

const ElementVectorDefinition* ElementVectorBase::GetDefinition() const {
  return def_.get();
}

ElementVector::ElementVector(const ElementVectorDefinition::CPtr& def)
    : ElementVectorBase(def) {
  Construct();
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

void ElementVector::Construct(){
  elements_.clear();
  for (int i = 0; i < def_->GetNumElements(); i++) {
    elements_.push_back(def_->GetElementDescription(i).MakeElement());
  }
}

ElementVectorWrapper::ElementVectorWrapper(
      const ElementVectorDefinition::CPtr& def,
      const ElementVectorDefinition::CPtr& in): ElementVectorBase(def), inDef_(in),
                                                element_vector_(nullptr),
                                                const_element_vector_(nullptr){
}

ElementVectorWrapper::~ElementVectorWrapper() {
}

ElementVectorWrapper& ElementVectorWrapper::operator=(const ElementVectorBase& other) {
  dynamic_cast<ElementVectorBase&>(*this) = other;
  return *this;
}

int ElementVectorWrapper::GetNumElement() const {
  return indexMap_.size();
}

void ElementVectorWrapper::ComputeMap() {
  indexMap_.resize(def_->GetNumElements());
  for (int i = 0; i < def_->GetNumElements(); i++) {
    indexMap_[i] = inDef_->FindName(def_->GetName(i));
    DLOG_IF(ERROR, indexMap_.at(i) == -1) << "Element name not found";
  }
}

void ElementVectorWrapper::SetElementVector(ElementVectorBase* element_vector) {
  element_vector->MatchesDefinition(*inDef_);
  element_vector_ = element_vector;
  const_element_vector_ = element_vector;
}

void ElementVectorWrapper::SetElementVector(const ElementVectorBase* element_vector) const {
  element_vector->MatchesDefinition(*inDef_);
  const_element_vector_ = element_vector;
}

void ElementVectorWrapper::EmbedJacobian(MatRef<> out,
                                        const MatCRef<>& in,
                                        int rowOffset) const {
  const int rows = in.rows();
  for (int i = 0; i < GetNumElement(); ++i) {
    const int cols = GetElement(i)->GetDim();
    out.block(rowOffset, const_element_vector_->GetStart(indexMap_.at(i)), rows, cols) =
        in.block(0, GetStart(i), rows, cols);
  }
}

}
