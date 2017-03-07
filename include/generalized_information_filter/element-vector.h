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
  bool MatchesDefinition(const ElementVectorDefinition& def) const;
  ElementVectorBase& operator=(const ElementVectorBase& other);
  virtual ElementBase* GetElement(int i) = 0;
  virtual const ElementBase* GetElement(int i) const = 0;
  template<typename T>
  inline T& GetValue(int i);
  template<typename T>
  inline const T& GetValue(int i) const;
  template<typename T>
  T& GetValue(const std::string& name);
  template<typename T>
  const T& GetValue(const std::string& name) const;
  inline int GetDim() const;
  virtual int GetNumElement() const = 0;
  inline int GetStart(int outer_index) const;
  inline int GetOuter(int i) const;
  inline int GetInner(int i) const;
  inline std::string GetName(int outer_index) const;
  inline int FindName(const std::string& name) const;
  std::string Print() const;
  void SetIdentity();
  void SetRandom();
  void BoxPlus(const VecCRefX& vec, ElementVectorBase* out) const;
  void BoxMinus(const ElementVectorBase& ref, VecRefX vec) const;
  const ElementVectorDefinition* GetDefinition() const;

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
  ElementVector& operator=(const ElementVector& other);
  ElementVector& operator=(ElementVector& other);
  int GetNumElement() const;
  inline ElementBase* GetElement(int i);
  inline const ElementBase* GetElement(int i) const;
  void Construct();

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
  inline ElementBase* GetElement(int i);
  inline const ElementBase* GetElement(int i) const;
  void ComputeMap();
  void SetElementVector(ElementVectorBase* element_vector);
  void SetElementVector(const ElementVectorBase* element_vector) const;
  void EmbedJacobian(MatRef<> out, const MatCRef<>& in, int rowOffset = 0) const;

 protected:
  ElementVectorBase* element_vector_;
  mutable const ElementVectorBase* const_element_vector_;
  const ElementVectorDefinition::CPtr inDef_;
  std::vector<int> indexMap_;
};

// ==================== Implementation ==================== //
template<typename T>
T& ElementVectorBase::GetValue(int i) {
  return GetElement(i)->template GetValue<T>();
}

template<typename T>
const T& ElementVectorBase::GetValue(int i) const {
  return GetElement(i)->template GetValue<T>();
}

template<typename T>
T& ElementVectorBase::GetValue(const std::string& name) {
  int i = FindName(name);
  DLOG_IF(ERROR, i == -1) << "Element name not found";
  return GetValue<T>(i);
}

template<typename T>
const T& ElementVectorBase::GetValue(const std::string& name) const {
  int i = FindName(name);
  DLOG_IF(ERROR, i == -1) << "Element name not found";
  return GetValue<T>(i);
}

int ElementVectorBase::GetDim() const {
  return GetDefinition()->GetDim();
}

int ElementVectorBase::GetStart(int outer_index) const {
  return GetDefinition()->GetStart(outer_index);
}

int ElementVectorBase::GetOuter(int i) const {
  return GetDefinition()->GetOuter(i);
}

int ElementVectorBase::GetInner(int i) const {
  return GetDefinition()->GetInner(i);
}

std::string ElementVectorBase::GetName(int outer_index) const {
  return GetDefinition()->GetName(outer_index);
}

int ElementVectorBase::FindName(const std::string& name) const {
  return GetDefinition()->FindName(name);
}

ElementBase* ElementVector::GetElement(int i) {
  return elements_.at(i).get();
}

const ElementBase* ElementVector::GetElement(int i) const {
  return elements_.at(i).get();
}

ElementBase* ElementVectorWrapper::GetElement(int i) {
  DLOG_IF(FATAL, element_vector_ == nullptr) << "Element vector not set for wrapper";
  return element_vector_->GetElement(indexMap_.at(i));
}

const ElementBase* ElementVectorWrapper::GetElement(int i) const {
  DLOG_IF(FATAL, const_element_vector_ == nullptr) << "Element vector not set for wrapper";
  return const_element_vector_->GetElement(indexMap_.at(i));
}

}
#endif /* GIF_ELEMENTVECTOR_HPP_ */
