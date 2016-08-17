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
  StateBase(const std::shared_ptr<const StateDefinition>& def): def_(def){};
  virtual ~StateBase(){};
  StateBase& operator=(const StateBase& other){
    for(int i=0;i<getNumElement();i++){
      *getElement(i) = *other.getElement(i);
    }
    return *this;
  }
  virtual std::shared_ptr<ElementBase> getElement(int i) = 0;
  virtual const std::shared_ptr<const ElementBase> getElement(int i) const = 0;
  template<typename T>
  inline T& getValue(int i){
    return std::dynamic_pointer_cast<Element<T>>(getElement(i))->get();
  }
  template<typename T>
  inline T& getValue(int i) const{
    return std::dynamic_pointer_cast<const Element<T>>(getElement(i))->get();
  }
  template<typename T>
  T& getValue(const std::string& name){
    int i = def_->findName(name);
    assert(i != -1);
    return getValue<T>(i);
  }
  template<typename T>
  T& getValue(const std::string& name) const{
    int i = def_->findName(name);
    assert(i != -1);
    return getValue<T>(i);
  }
  inline int getDim() const{
    return def_->getDim();
  }
  inline int getNumElement() const{
    return def_->getNumElement();
  }
  inline int getStart(int i) const{
    return def_->getStart(i);
  }
  inline int getOuter(int i) const{
    return def_->getOuter(i);
  }
  inline int getInner(int i) const{
    return def_->getInner(i);
  }
  void print() const{
    for(int i=0;i<getNumElement();i++){
      getElement(i)->print();
    }
  }
  void setIdentity(){
    for(int i=0;i<getNumElement();i++){
      getElement(i)->setIdentity();
    }
  }
  void boxplus(const VXD& vec, const std::shared_ptr<StateBase>& out) const{
    for(int i=0;i<getNumElement();i++){
      getElement(i)->boxplus(vec.block(getStart(i),0,getElement(i)->getDim(),1),out->getElement(i));
    }
  }
  void boxminus(const std::shared_ptr<const StateBase>& ref, VXD& vec) const{
    for(int i=0;i<getNumElement();i++){
      getElement(i)->boxminus(ref->getElement(i),vec.block(getStart(i),0,getElement(i)->getDim(),1));
    }
  }

 protected:
  const std::shared_ptr<const StateDefinition> def_;
};

class State: public StateBase{
 public:
  State(const std::shared_ptr<const StateDefinition>& def): StateBase(def){
    for(int i=0;i<def_->getNumElement();i++){
      elements_.push_back(def_->getElementDefinition(i)->newElement());
    }
  };
  virtual ~State(){};
  State& operator=(const StateBase& other){
    dynamic_cast<const StateBase&>(*this) = other;
    return *this;
  }
  inline std::shared_ptr<ElementBase> getElement(int i){
    return elements_.at(i);
  }
  inline const std::shared_ptr<const ElementBase> getElement(int i) const{
    return elements_.at(i);
  }

 protected:
  std::vector<std::shared_ptr<ElementBase>> elements_;
};

class StateWrapper: public StateBase{ // TODO: implement definiton checking (Wrapper-State, Wrapper-Definition)
 public:
  StateWrapper(const std::shared_ptr<const StateDefinition>& def, const std::shared_ptr<const StateDefinition>& in): StateBase(def), in_(in){
    computeMap();
  };
  ~StateWrapper(){};
  StateWrapper& operator=(const StateBase& other){
    dynamic_cast<const StateBase&>(*this) = other;
    return *this;
  }
  inline std::shared_ptr<ElementBase> getElement(int i){
    return state_->getElement(indexMap_[i]);
  }
  inline const std::shared_ptr<const ElementBase> getElement(int i) const{
    return constState_->getElement(indexMap_[i]);
  }
  void computeMap(){
    indexMap_.resize(def_->getNumElement());
    for(int i=0;i<def_->getNumElement();i++){
      indexMap_[i] = in_->findName(def_->getName(i));
      assert(indexMap_[i] != -1);
    }
  }
  void setState(const std::shared_ptr<StateBase>& state){
    state_ = state;
  }
  void setState(const std::shared_ptr<const StateBase>& state) const{
    constState_ = state;
  }

 protected:
  std::shared_ptr<StateBase> state_;
  mutable std::shared_ptr<const StateBase> constState_;
  const std::shared_ptr<const StateDefinition> in_;
  std::vector<int> indexMap_;
};

}

#endif /* GIF_STATE_HPP_ */
