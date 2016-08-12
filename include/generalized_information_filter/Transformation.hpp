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

template<typename PackOut, typename PackIn>
class Transformation;

template<typename... Out, typename... In>
class Transformation<ElementPack<Out...>,ElementPack<In...>>: public Model<Transformation<ElementPack<Out...>,ElementPack<In...>>, ElementPack<Out...>, ElementPack<In...>>{
 public:
  using mtBase = Model<Transformation<ElementPack<Out...>,ElementPack<In...>>, ElementPack<Out...>, ElementPack<In...>>;
  typedef Transformation<ElementPack<Out...>,ElementPack<In...>> mtTransformation;
  Transformation(const std::array<std::string,mtBase::n_>& namesOut, const std::array<std::string,mtBase::m_>& namesIn):
      mtBase(namesOut, std::forward_as_tuple(namesIn)), J_(0,0){};
  virtual ~Transformation(){};
  virtual void evalTransform(Out&... outs, const In&... ins) = 0;
  virtual void jacTransform(MXD& J, const In&... ins) = 0;
  void eval(Out&... outs, const In&... ins){
    evalTransform(outs...,ins...);
  }
  template<int j, typename std::enable_if<(j==0)>::type* = nullptr>
  void jac(MXD& J, const In&... ins){
    jacTransform(J,ins...);
  }
  void transformState(Out&... outs, const In&... ins){ // TODO: switch to State*
    evalTransform(outs..., ins...);
  }
  void transformCovMat(const In&... ins,const MXD& inputCov, MXD& outputCov){
    jacTransform(J_,ins...);
    outputCov = J_*inputCov*J_.transpose();
    postProcess(outputCov, ins...);
  }
  virtual void postProcess(MXD& P, const In&... ins){};
  std::shared_ptr<StateDefinition> inputDefinition(){
    return this->inDefinitions_[0];
  }
  std::shared_ptr<StateDefinition> outputDefinition(){
    return this->outDefinition_;
  }

 protected:
  MXD J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
