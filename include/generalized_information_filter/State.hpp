/*
 * State.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_STATE_HPP_
#define GIF_STATE_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Element.hpp"
#include "generalized_information_filter/StateDefinition.hpp"

namespace GIF{

class StateBase{
 public:
  StateBase(const std::shared_ptr<const StateDefinition>& def);
  virtual ~StateBase();
  bool matchesDef(const std::shared_ptr<const StateDefinition>& def) const;
  StateBase& operator=(const StateBase& other);
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
  void boxplus(const VXD& vec, const std::shared_ptr<StateBase>& out) const;
  void boxminus(const std::shared_ptr<const StateBase>& ref, VXD& vec) const;

 protected:
  const std::shared_ptr<const StateDefinition> def_;
};

class State: public StateBase{
 public:
  State(const std::shared_ptr<const StateDefinition>& def);
  virtual ~State();
  State& operator=(const StateBase& other);
  int getNumElement() const;
  inline std::shared_ptr<ElementBase> getElement(int i);
  inline std::shared_ptr<const ElementBase> getElement(int i) const;

 protected:
  std::vector<std::shared_ptr<ElementBase>> elements_;
};

class StateWrapper: public StateBase{ // TODO: implement definiton checking (Wrapper-State, Wrapper-Definition)
 public:
  StateWrapper(const std::shared_ptr<const StateDefinition>& def, const std::shared_ptr<const StateDefinition>& in);
  ~StateWrapper();
  StateWrapper& operator=(const StateBase& other);
  int getNumElement() const;
  inline std::shared_ptr<ElementBase> getElement(int i);
  inline std::shared_ptr<const ElementBase> getElement(int i) const;
  void computeMap();
  void setState(const std::shared_ptr<StateBase>& state);
  void setState(const std::shared_ptr<const StateBase>& state) const;

 protected:
  std::shared_ptr<StateBase> state_;
  mutable std::shared_ptr<const StateBase> constState_;
  const std::shared_ptr<const StateDefinition> in_;
  std::vector<int> indexMap_;
};


// ==================== Implementation ==================== //
template<typename T>
inline T& StateBase::getValue(int i){
  return std::dynamic_pointer_cast<Element<T>>(getElement(i))->get();
}

template<typename T>
inline T& StateBase::getValue(int i) const{
  return std::dynamic_pointer_cast<const Element<T>>(getElement(i))->get();
}

template<typename T>
T& StateBase::getValue(const std::string& name){
  assert(matchesDef(def_));
  int i = def_->findName(name);
  assert(i != -1);
  return getValue<T>(i);
}

template<typename T>
T& StateBase::getValue(const std::string& name) const{
  assert(matchesDef(def_));
  int i = def_->findName(name);
  assert(i != -1);
  return getValue<T>(i);
}

int StateBase::getDim() const{
  assert(matchesDef(def_));
  return def_->getDim();
}

int StateBase::getStart(int i) const{
  assert(matchesDef(def_));
  return def_->getStart(i);
}

int StateBase::getOuter(int i) const{
  assert(matchesDef(def_));
  return def_->getOuter(i);
}

int StateBase::getInner(int i) const{
  assert(matchesDef(def_));
  return def_->getInner(i);
}

std::shared_ptr<ElementBase> State::getElement(int i){
  return elements_.at(i);
}

std::shared_ptr<const ElementBase> State::getElement(int i) const{
  return elements_.at(i);
}

std::shared_ptr<ElementBase> StateWrapper::getElement(int i){
  return state_->getElement(indexMap_[i]);
}

std::shared_ptr<const ElementBase> StateWrapper::getElement(int i) const{
  return constState_->getElement(indexMap_[i]);
}

}

#endif /* GIF_STATE_HPP_ */
