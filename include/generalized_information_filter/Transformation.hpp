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
  void transformState(const std::vector<const State*>& in, State* out){
    eval(in, out);
  }
  void transformCovMat(const std::vector<const State*>& in,const MXD& inputCov, MXD& outputCov){
    jac(in,J_);
    outputCov = J_*inputCov*J_.transpose();
    postProcess(in,outputCov);
  }
  virtual void postProcess(const std::vector<const State*>& in, MXD& cov){};

 private:
  MXD J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
