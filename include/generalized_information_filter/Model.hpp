/*
 * Model.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include <list>
#include <initializer_list>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/StateDefinition.hpp"

namespace GIF{

class Model{
 public:
  Model(){
    outputDefinition_ = nullptr;
  };
  Model(StateDefinition* outputDefinition, std::initializer_list<StateDefinition*> inputList): outputDefinition_(outputDefinition), inputDefinition_(inputList){};
  virtual ~Model(){};
  virtual void eval(const std::vector<const State*>& in, State* out) = 0;
  virtual void jac(const std::vector<const State*>& in, MXD& out) = 0;
  void jacFD(const std::vector<const State*>& in, MXD& J, const double& delta, const int& c = 0){
    J.setZero();
    State* stateDis = inputDefinition_[c]->newState();
    *stateDis = *in[c];
    std::vector<const State*> inDis(in);
    inDis[c] = stateDis;
    State* outRef = outputDefinition_->newState();
    State* outDis = outputDefinition_->newState();
    eval(inDis,outRef);
    VXD difIn(inputDefinition_[c]->getDim());
    VXD difOut(outputDefinition_->getDim());
    for(int i=0; i<inputDefinition_[c]->getDim(); i++){
      difIn.setZero();
      difIn(i) = delta;
      inputDefinition_[c]->boxplus(in[c],difIn,stateDis);
      eval(inDis,outDis);
      outputDefinition_->boxminus(outDis,outRef,difOut);
      J.col(i) = difOut/delta;
    }

    delete stateDis;
    delete outRef;
    delete outDis;
  }

 protected:
  std::vector<StateDefinition*> inputDefinition_;
  StateDefinition* outputDefinition_;
};

}

#endif /* GIF_MODEL_HPP_ */
