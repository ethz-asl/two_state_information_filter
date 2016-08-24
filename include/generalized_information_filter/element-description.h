#ifndef GIF_ELEMENTDESCRIPTION_HPP_
#define GIF_ELEMENTDESCRIPTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element.h"

namespace GIF {

class ElementDescriptionBase {
 public:
  ElementDescriptionBase() {}

  virtual ~ElementDescriptionBase() {}

  virtual std::shared_ptr<ElementBase> newElement() const = 0;
  virtual bool isSameDef(
      const std::shared_ptr<const ElementDescriptionBase>& in) const = 0;
  virtual bool isOfDef(const std::shared_ptr<const ElementBase>& in) const = 0;
  virtual int getDim() const = 0;
};

template<typename T>
class ElementDescription : public ElementDescriptionBase {
 public:
  ElementDescription() {
  }
  ~ElementDescription() {
  }
  std::shared_ptr<ElementBase> newElement() const {
    return std::shared_ptr < Element < T >> (new Element<T>(this));
  }
  bool isSameDef(const std::shared_ptr<const ElementDescriptionBase>& in) const {
    return std::dynamic_pointer_cast<const ElementDescription<T>>(in).get()
        != nullptr;
  }
  bool isOfDef(const std::shared_ptr<const ElementBase>& in) const {
    return std::dynamic_pointer_cast<const Element<T>>(in).get() != nullptr;
  }
  inline int getDim() const {
    return ElementTraits<T>::d_;
  }
};

}

#endif /* GIF_ELEMENTDESCRIPTION_HPP_ */
