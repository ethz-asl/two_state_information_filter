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

class State{
 public:
  State(const std::shared_ptr<const StateDefinition> def): def_(def){
    for(int i=0;i<def_->getNumElement();i++){
      elements_.push_back(def_->getElementDefinition(i)->newElement());
    }
  };
  virtual ~State(){};
  State& operator=(const State& other){
    for(int i=0;i<elements_.size();i++){
      *getElement(i) = *other.getElement(i);
    }
    return *this;
  }
  inline std::shared_ptr<ElementBase>& getElement(int i){
    return elements_.at(i);
  }
  inline const std::shared_ptr<const ElementBase> getElement(int i) const{
    return elements_.at(i);
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
    for(auto e : elements_){
      e->print();
    }
  }
  void init(){
    for(auto e : elements_){
      e->setIdentity();
    }
  }
  void boxplus(const VXD& vec, const std::shared_ptr<State> out) const{
    for(int i=0;i<getNumElement();i++){
      getElement(i)->boxplus(vec.block(getStart(i),0,getElement(i)->getDim(),1),out->getElement(i));
    }
  }
  void boxminus(const std::shared_ptr<const State> ref, VXD& vec) const{
    for(int i=0;i<getNumElement();i++){
      getElement(i)->boxminus(ref->getElement(i),vec.block(getStart(i),0,getElement(i)->getDim(),1));
    }
  }

 protected:
  std::vector<std::shared_ptr<ElementBase>> elements_;
  const std::shared_ptr<const StateDefinition> def_;
};

//class StateWrapper{ // TODO: implement definiton checking (Wrapper-State, Wrapper-Definition)
// public:
//  StateWrapper(){};
//  ~StateWrapper(){};
//  void computeMap(){
//    indexMap_.resize(outDefinition_->namesMap_.size());
//    for(auto outElement : outDefinition_->namesMap_){
//      indexMap_.at(outElement.second) = inDefinition_->namesMap_.at(outElement.first);
//    }
//  }
//  void wrap(const std::shared_ptr<State>& out, const std::shared_ptr<State>& in){
//    for(int i=0;i<indexMap_.size();i++){
//      out->getElement(i) = in->getElement(indexMap_[i]);
//    }
//  }
//
// protected:
//  std::shared_ptr<StateDefinition> inDefinition_;
//  std::shared_ptr<StateDefinition> outDefinition_;
//  std::vector<int> indexMap_;
//};

}

#endif /* GIF_STATE_HPP_ */
