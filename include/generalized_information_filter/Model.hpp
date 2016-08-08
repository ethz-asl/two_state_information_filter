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
  Model(){};
  virtual ~Model(){};
  virtual void _eval(const std::vector<const State*>& in, State* out) = 0;
  virtual void _jac(const std::vector<const State*>& in, MXD& out, int c) = 0;
  void _jacFD(const std::vector<const State*>& in, MXD& J, int c, const double& delta){
    J.setZero();
    State* stateDis = inputDefinitions_[c]->newState();
    *stateDis = *in[c];
    std::vector<const State*> inDis(in);
    inDis[c] = stateDis;
    State* outRef = outputDefinition_->newState();
    State* outDis = outputDefinition_->newState();
    _eval(inDis,outRef);
    VXD difIn(inputDefinitions_[c]->getDim());
    VXD difOut(outputDefinition_->getDim());
    for(int i=0; i<inputDefinitions_[c]->getDim(); i++){
      difIn.setZero();
      difIn(i) = delta;
      inputDefinitions_[c]->boxplus(in[c],difIn,stateDis);
      _eval(inDis,outDis);
      outputDefinition_->boxminus(outDis,outRef,difOut);
      J.col(i) = difOut/delta;
    }
    delete stateDis;
    delete outRef;
    delete outDis;
  }
  void initStateDefinitions(std::shared_ptr<StateDefinition> outputDefinition, std::initializer_list<std::shared_ptr<StateDefinition>> inputList){
    outputDefinition_ = outputDefinition;
    inputDefinitions_ = inputList;
    buildStateDefinitions();
  }

 protected:
  virtual void buildStateDefinitions() = 0;

 private:
  std::vector<std::shared_ptr<StateDefinition>> inputDefinitions_;
  std::shared_ptr<StateDefinition> outputDefinition_;
};

}

#endif /* GIF_MODEL_HPP_ */
