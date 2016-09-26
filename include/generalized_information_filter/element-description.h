#ifndef GIF_ELEMENTDESCRIPTION_HPP_
#define GIF_ELEMENTDESCRIPTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"

namespace GIF {

/*! \brief Element Description Base
 *         Base class for element descriptions. Used to  implement the element
 *         vector definition (vector of element descriptions).
 */
class ElementDescriptionBase {
 public:
  typedef std::shared_ptr<ElementDescriptionBase> Ptr;
  typedef std::shared_ptr<const ElementDescriptionBase> CPtr;
  ElementDescriptionBase() {}
  virtual ~ElementDescriptionBase() {}

  virtual ElementBase::Ptr MakeElement(const CPtr& description) const = 0;
  virtual bool MatchesDescription(const ElementDescriptionBase::CPtr& in) const = 0;
  virtual bool MatchesDescription(const ElementBase& in) const = 0;
  virtual int GetDim() const = 0;
  virtual bool IsVectorSpace() const = 0;
};

/*! \brief Templated form of element descriptions.
 *         Implements the virtual methods of the base class (mainly based on
 *         dynamic casting).
 */
template<typename T>
class ElementDescription : public ElementDescriptionBase {
 public:
  typedef std::shared_ptr<ElementDescription<T>> Ptr;
  typedef std::shared_ptr<const ElementDescription<T>> CPtr;
  ElementDescription() {
  }
  ~ElementDescription() {
  }
  ElementBase::Ptr MakeElement(const ElementDescriptionBase::CPtr& description) const {
    return std::make_shared<Element<T>>(
        std::dynamic_pointer_cast<const ElementDescription<T>>(description));
  }
  bool MatchesDescription(const ElementDescriptionBase::CPtr& in) const {
    return (bool)std::dynamic_pointer_cast<const ElementDescription<T>>(in);
  }
  bool MatchesDescription(const ElementBase& in) const {
    const ElementBase* inPtr = &in;
    return dynamic_cast<const Element<T>*>(inPtr) != nullptr;
  }
  inline int GetDim() const {
    return ElementTraits<T>::kDim;
  }
  inline bool IsVectorSpace() const {
    return ElementTraits<T>::kIsVectorSpace;
  }
};

}

#endif /* GIF_ELEMENTDESCRIPTION_HPP_ */
