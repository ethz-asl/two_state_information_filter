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

  virtual ElementBase::Ptr MakeElement() const = 0;
  virtual bool MatchesDescription(const ElementDescriptionBase& in) const = 0;
  virtual bool MatchesDescription(const ElementBase& in) const = 0;
  virtual int GetDimension() const = 0;
  virtual ElementDescriptionBase::Ptr Copy() const = 0;
};

/*! \brief Templated form of element descriptions.
 *         Implements the virtual methods of the base class (mainly based on
 *         dynamic casting).
 */
template<typename T>
class ElementDescription : public ElementDescriptionBase {
 public:
  ElementDescription() {
  }
  ~ElementDescription() {
  }
  ElementBase::Ptr MakeElement() const {
    return std::make_shared<Element<T>>(this);
  }
  bool MatchesDescription(const ElementDescriptionBase& in) const {
    const ElementDescriptionBase* inPtr = &in;
    return dynamic_cast<const ElementDescription<T>*>(inPtr) != nullptr;
  }
  bool MatchesDescription(const ElementBase& in) const {
    const ElementBase* inPtr = &in;
    return dynamic_cast<const Element<T>*>(inPtr) != nullptr;
  }
  inline int GetDimension() const {
    return ElementTraits<T>::d_;
  }
  ElementDescriptionBase::Ptr Copy() const{
    return std::make_shared<ElementDescription<T>>();
  }
};

}

#endif /* GIF_ELEMENTDESCRIPTION_HPP_ */
