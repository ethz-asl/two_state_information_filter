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
  Transformation(StateDefinition* inputDefinition, StateDefinition* outputDefinition): Model(outputDefinition,{inputDefinition}){};
  virtual ~Transformation(){};
  virtual void eval(const State* in, State* out) = 0;
  virtual void jac(const State* in, MXD& out) = 0;
  virtual void _eval(const std::vector<const State*>& in, State* out){
    eval(in[0],out);
  }
  virtual void _jac(const std::vector<const State*>& in, MXD& out){
    jac(in[0],out);
  }
  void jacFD(const State* in, MXD& out, const double& delta){
    const std::vector<const State*> inVec({in});
    _jacFD(inVec,out,delta,0);
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

 private:
  MXD J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
