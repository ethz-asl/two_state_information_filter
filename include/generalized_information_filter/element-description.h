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
  virtual bool MatchesDescription(
      const ElementDescriptionBase::CPtr& in) const = 0;
  virtual bool MatchesDescription(const ElementBase::CPtr& in) const = 0;
  virtual int GetDimension() const = 0;
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
    return std::shared_ptr<Element<T>>(new Element<T>(this));
  }
  bool MatchesDescription(const ElementDescriptionBase::CPtr& in) const {
    return std::dynamic_pointer_cast<const ElementDescription<T>>(in).get()
        != nullptr;
  }
  bool MatchesDescription(const ElementBase::CPtr& in) const {
    return std::dynamic_pointer_cast<const Element<T>>(in).get() != nullptr;
  }
  inline int GetDimension() const {
    return ElementTraits<T>::d_;
  }
};

}

#endif /* GIF_ELEMENTDESCRIPTION_HPP_ */
