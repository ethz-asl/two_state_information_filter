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
      mtBase(namesOut, std::forward_as_tuple(namesIn)), J_((int)ElementPack<Out...>::d_,(int)ElementPack<In...>::d_){};
  virtual ~Transformation(){};

  // User implementations
  virtual void evalTransform(Out&... outs, const In&... ins) const = 0;
  virtual void jacTransform(MXD& J, const In&... ins) const  = 0;

  // Wrapping from user interface to base
  void jacFD(MXD& J, const std::shared_ptr<const StateBase>& in, const double& delta = 1e-6){
    const std::array<std::shared_ptr<const StateBase>,1> ins = {in};
    this->template _jacFD<0>(J,ins,delta);
  }
  void transformState(const std::shared_ptr<StateBase>& out, const std::shared_ptr<const StateBase>& in){
    const std::array<std::shared_ptr<const StateBase>,1> ins = {in};
    this->template _eval(out,ins);
  }
  void transformCovMat(MXD& outputCov, const std::shared_ptr<const StateBase>& in,const MXD& inputCov){
    const std::array<std::shared_ptr<const StateBase>,1> ins = {in};
    this->template _jac<0>(J_,ins);
    outputCov = J_*inputCov*J_.transpose();
  }
  template<int n, int m>
  void setJacBlock(MXD& J, const Eigen::Matrix<double,ElementPack<Out...>::template getDim<n>(),ElementPack<In...>::template getDim<m>()>& B) const{
    this->template _setJacBlock<0,n,m>(J,B);
  }
  bool testJac(const std::shared_ptr<const StateBase>& in, const double& delta = 1e-6, const double& th = 1e-6) const{
    const std::array<std::shared_ptr<const StateBase>,1> ins = {in};
    return this->template _testJacInput<0>(ins,delta,th);
  }

  // Access to definitions
  std::shared_ptr<StateDefinition> outputDefinition(){
    return this->outDefinition_;
  }
  std::shared_ptr<StateDefinition> inputDefinition(){
    return this->inDefinitions_[0];
  }

 protected:
  friend mtBase;

  // Wrapping from base to user implementation
  void eval(Out&... outs, const In&... ins) const{
    evalTransform(outs...,ins...);
  }
  template<int j, typename std::enable_if<(j==0)>::type* = nullptr>
  void jac(MXD& J, const In&... ins) const{
    jacTransform(J,ins...);
  }

  MXD J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
