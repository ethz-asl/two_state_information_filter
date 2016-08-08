/*
 * Transformation.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_TRANSFORMATION_HPP_
#define GIF_TRANSFORMATION_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Model.hpp"

namespace GIF{

class Transformation: public Model{
 public:
  Transformation(): inputDefinition_(new StateDefinition()), outputDefinition_(new StateDefinition()), J_(0,0){};
  virtual ~Transformation(){};
  virtual void eval(const State* in, State* out) = 0;
  virtual void jac(const State* in, MXD& out) = 0;
  virtual void _eval(const std::vector<const State*>& in, State* out){
    eval(in[0],out);
  }
  virtual void _jac(const std::vector<const State*>& in, MXD& out, int c){
    assert(c==0);
    jac(in[0],out);
  }
  void jacFD(const State* in, MXD& out, const double& delta){
    const std::vector<const State*> inVec({in});
    _jacFD(inVec,out,0,delta);
  }
  void transformState(const State* in, State* out){
    eval(in, out);
  }
  void transformCovMat(const State* in,const MXD& inputCov, MXD& outputCov){
    jac(in,J_);
    outputCov = J_*inputCov*J_.transpose();
    postProcess(in,outputCov);
  }
  virtual void postProcess(const State* in, MXD& cov){};
  void initStateDefinitions(){
    Model::initStateDefinitions(outputDefinition_,{inputDefinition_});
  }
  std::shared_ptr<const StateDefinition> inputDefinition(){
    return inputDefinition_;
  }
  std::shared_ptr<const StateDefinition> outputDefinition(){
    return outputDefinition_;
  }

 protected:
  std::shared_ptr<StateDefinition> inputDefinition_;
  std::shared_ptr<StateDefinition> outputDefinition_;
  MXD J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
