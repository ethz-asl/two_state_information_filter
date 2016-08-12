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

namespace GIF{

class State{
 public:
  State(){
    d_ = 0;
  };
  virtual ~State(){};
  State& operator=(const State& other){
    for(int i=0;i<elements_.size();i++){
      *elements_.at(i).first = *other.elements_.at(i).first;
    }
    return *this;
  }
  void addElement(const std::shared_ptr<ElementBase>& e){
    elements_.push_back(std::pair<std::shared_ptr<ElementBase>,int>(e,d_));
    d_ += e->getDim();
  }
  std::shared_ptr<ElementBase>& getElement(int i){
    return std::get<0>(elements_.at(i));
  }
  const std::shared_ptr<const ElementBase> getElement(int i) const{
    return std::get<0>(elements_.at(i));
  }
  const int& getIndex(int i) const{
    return std::get<1>(elements_.at(i));
  }
  int getDim() const{
    return d_;
  }
  void print() const{
    for(auto e : elements_){
      e.first->print();
    }
  }
  void init(){
    for(auto e : elements_){
      e.first->init();
    }
  }
  void boxplus(const VXD& vec, const std::shared_ptr<State> out) const{
    for(int i=0;i<elements_.size();i++){
      getElement(i)->boxplus(vec.block(getIndex(i),0,getElement(i)->getDim(),1),out->getElement(i));
    }
  }
  void boxminus(const std::shared_ptr<const State> ref, VXD& vec) const{
    for(int i=0;i<elements_.size();i++){
      getElement(i)->boxminus(ref->getElement(i),vec.block(getIndex(i),0,getElement(i)->getDim(),1));
    }
  }
  int getOuter(int ind) const{
    int i=0;
    while(ind >= getElement(i)->getDim()){
      ind -= getElement(i)->getDim();
      ++i;
    }
    return i;
  }
  int getInner(int ind) const{
    return ind-getIndex(getOuter(ind));
  }

 protected:
  std::vector<std::pair<std::shared_ptr<ElementBase>,int>> elements_;
  int d_;
};

}

#endif /* GIF_STATE_HPP_ */
