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
  }
  State(const State& other){
    d_ = other.d_;;
    for(auto e : other.elements_){
      elements_.push_back(e->clone());
    }
  }
  virtual ~State(){
    for(auto e : elements_){
      delete e;
    }
  };
  virtual State* clone() const = 0;
  State& operator=(const State& other){
    assert(this->d_==other.d_);
    for(int i=0;i<elements_.size();i++){
      *elements_[i] = *other.elements_[i];
    }
    return *this;
  }
  template<typename T>
  Element<T>* addElement(){
    Element<T>* e = new Element<T>();
    e->init();
    e->i_ = d_;
    d_ += e->d_;
    elements_.push_back(e);
    return e;
  }
  ElementBase& getElement(int ID){
    return *elements_[ID];
  }
  const ElementBase& getElement(int ID) const{
    return *elements_[ID];
  }
  void print() const{
    for(auto e : elements_){
      e->print();
    }
  }
  void init() const{
    for(auto e : elements_){
      e->init();
    }
  }
  void boxplus(const VXD& vec, State& out) const{
    int c=0;
    for(int i=0;i < elements_.size();i++){
      elements_[i]->boxplus(vec.block(c,0,elements_[i]->d_,1),out.elements_[i]);
      c+=elements_[i]->d_;
    }
  }
  void boxminus(const State& ref, VXD& vec) const{
    int c=0;
    for(int i=0;i < elements_.size();i++){
      elements_[i]->boxminus(ref.elements_[i],vec.block(c,0,elements_[i]->d_,1));
      c+=elements_[i]->d_;
    }
  }
  int d_;
 private:
  std::vector<ElementBase*> elements_;
};

}

#endif /* GIF_STATE_HPP_ */
