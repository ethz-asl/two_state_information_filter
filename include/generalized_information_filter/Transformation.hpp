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
  Transformation(): outputDefinition_(new StateDefinition()), inputDefinition_(new StateDefinition()), J_(0,0){
    ElementPack<Out...>::addElementToDefinition(outputDefinition_);
    ElementPack<In...>::addElementToDefinition(inputDefinition_);
  };
  virtual ~Transformation(){};
  virtual void eval(Out&... outs, const In&... ins) = 0;
  virtual void jac(MXD* J, const In&... ins) = 0;
  void transformState(Out&... outs, const In&... ins){
    eval(outs..., ins...);
  }
  void transformCovMat(const In&... ins,const MXD& inputCov, MXD& outputCov){
    jac(J_,ins...);
    outputCov = J_*inputCov*J_.transpose();
    postProcess(ins...,outputCov);
  }
  virtual void postProcess(const In&... ins, MXD& cov){};
  std::shared_ptr<StateDefinition> inputDefinition(){
    return inputDefinition_;
  }
  std::shared_ptr<StateDefinition> outputDefinition(){
    return outputDefinition_;
  }

 protected:
  std::shared_ptr<StateDefinition> inputDefinition_;
  std::shared_ptr<StateDefinition> outputDefinition_;
  MXD J_;
};

//class Transformation: public Model{
// public:
//  Transformation(): inputDefinition_(new StateDefinition()), outputDefinition_(new StateDefinition()), J_(0,0){};
//  virtual ~Transformation(){};
//  virtual void eval(const State* in, State* out) = 0;
//  virtual void jac(const State* in, MXD& out) = 0;
//  virtual void _eval(const std::vector<const State*>& in, State* out){
//    eval(in[0],out);
//  }
//  virtual void _jac(const std::vector<const State*>& in, MXD& out, int c){
//    assert(c==0);
//    jac(in[0],out);
//  }
//  void jacFD(const State* in, MXD& out, const double& delta){
//    const std::vector<const State*> inVec({in});
//    _jacFD(inVec,out,0,delta);
//  }
//  void transformState(const State* in, State* out){
//    eval(in, out);
//  }
//  void transformCovMat(const State* in,const MXD& inputCov, MXD& outputCov){
//    jac(in,J_);
//    outputCov = J_*inputCov*J_.transpose();
//    postProcess(in,outputCov);
//  }
//  virtual void postProcess(const State* in, MXD& cov){};
//  void initStateDefinitions(){
//    Model::initStateDefinitions(outputDefinition_,{inputDefinition_});
//  }
//  std::shared_ptr<const StateDefinition> inputDefinition(){
//    return inputDefinition_;
//  }
//  std::shared_ptr<const StateDefinition> outputDefinition(){
//    return outputDefinition_;
//  }
//
// protected:
//  std::shared_ptr<StateDefinition> inputDefinition_;
//  std::shared_ptr<StateDefinition> outputDefinition_;
//  MXD J_;
//};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
