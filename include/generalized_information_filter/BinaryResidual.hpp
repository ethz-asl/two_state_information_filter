/*
 * BinaryResidual.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_BINARYRESIDUAL_HPP_
#define GIF_BINARYRESIDUAL_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Model.hpp"
#include "generalized_information_filter/Measurement.hpp"

namespace GIF{

class BinaryResidualBase{
 public:
  BinaryResidualBase(bool isUnary = false, bool isSplitable = false, bool isMergeable = false): isUnary_(isUnary), isSplitable_(isSplitable), isMergeable_(isMergeable){
  };
  virtual ~BinaryResidualBase(){};
  virtual void evalResidual(const std::shared_ptr<StateBase>& res, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const = 0;
  virtual void jacPre(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const  = 0;
  virtual void jacPos(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const  = 0;
  virtual void jacNoi(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const  = 0;
  virtual std::shared_ptr<StateDefinition> resDefinition() const = 0;
  virtual std::shared_ptr<StateDefinition> preDefinition() const = 0;
  virtual std::shared_ptr<StateDefinition> posDefinition() const = 0;
  virtual std::shared_ptr<StateDefinition> noiDefinition() const = 0;
  virtual void splitMeasurements(const std::shared_ptr<const MeasurementBase>& in, const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                                std::shared_ptr<const MeasurementBase>& out1, std::shared_ptr<const MeasurementBase>& out2) const = 0;
  virtual void mergeMeasurements(const std::shared_ptr<const MeasurementBase>& in1, const std::shared_ptr<const MeasurementBase>& in2,
                                const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, std::shared_ptr<const MeasurementBase>& out) const = 0;
  virtual void setMeas(const std::shared_ptr<const MeasurementBase>& meas) = 0;

  const bool isUnary_;
  const bool isSplitable_;
  const bool isMergeable_;
};

template<typename PackRes, typename PackPre, typename PackPos, typename PackNoi, typename Meas>
class BinaryResidual;

template<typename... Res, typename... Pre, typename... Pos, typename... Noi, typename Meas>
class BinaryResidual<ElementPack<Res...>,ElementPack<Pre...>,ElementPack<Pos...>,ElementPack<Noi...>,Meas>:
    public Model<BinaryResidual<ElementPack<Res...>,ElementPack<Pre...>,ElementPack<Pos...>,ElementPack<Noi...>,Meas>,
                                ElementPack<Res...>,ElementPack<Pre...>,ElementPack<Pos...>,ElementPack<Noi...>>,
    public BinaryResidualBase{
 public:
  typedef BinaryResidual<ElementPack<Res...>,ElementPack<Pre...>,ElementPack<Pos...>,ElementPack<Noi...>,Meas> mtBinaryRedidual;
  using mtBase = Model<mtBinaryRedidual,ElementPack<Res...>,ElementPack<Pre...>,ElementPack<Pos...>,ElementPack<Noi...>>;
  BinaryResidual(const std::array<std::string,ElementPack<Res...>::n_>& namesRes,
                 const std::array<std::string,ElementPack<Pre...>::n_>& namesPre,
                 const std::array<std::string,ElementPack<Pos...>::n_>& namesPos,
                 const std::array<std::string,ElementPack<Noi...>::n_>& namesNoi,
                 bool isUnary = false, bool isSplitable = false, bool isMergeable = false):
      mtBase(namesRes, std::forward_as_tuple(namesPre,namesPos,namesNoi)),
      BinaryResidualBase(isUnary,isSplitable,isMergeable){};
  virtual ~BinaryResidual(){};

  // Set measurement
  void setMeas(const std::shared_ptr<const MeasurementBase>& meas){
    meas_ = std::dynamic_pointer_cast<const Meas>(meas);
    if(!meas_){
      std::cout << "ERROR: Passing wrong measurement type" << std::endl;
    }
  }

  void splitMeasurements(const std::shared_ptr<const MeasurementBase>& in, const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                                std::shared_ptr<const MeasurementBase>& out1, std::shared_ptr<const MeasurementBase>& out2) const{
    // carefull: in/out1/out2 may point to the same elements
    if(isSplitable_){
      out1 = in;
      out2 = in;
    } else {
      std::cout << "ERROR: splitting of specific residual not supported/implemented!" << std::endl;
    }
  }

  void mergeMeasurements(const std::shared_ptr<const MeasurementBase>& in1, const std::shared_ptr<const MeasurementBase>& in2,
                                const TimePoint& t0, const TimePoint& t1, const TimePoint& t2, std::shared_ptr<const MeasurementBase>& out) const{
    if(isMergeable_){
      std::shared_ptr<MeasurementBase> newMeas(new Meas());
      VXD diff(in1->getDim());
      in1->boxminus(in2,diff);
      in2->boxplus(toSec(t1-t0)/toSec(t2-t0)*diff,newMeas);
      out = newMeas;
    } else {
      std::cout << "ERROR: merging of specific residual not supported!/implemented" << std::endl;
    }
  }

  // User implementations
  virtual void evalResidual(Res&... res, const Pre&... pre, const Pos&... pos, const Noi&... noi) const = 0;
  virtual void jacPre(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const  = 0;
  virtual void jacPos(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const  = 0;
  virtual void jacNoi(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const  = 0;

  // Wrapping from user interface to base
  void evalResidual(const std::shared_ptr<StateBase>& res, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const{
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->_eval(res,ins);
  }
  void jacPre(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const{
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jac<0>(J,ins);
  }
  void jacPos(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const{
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jac<1>(J,ins);
  }
  void jacNoi(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi) const{
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jac<2>(J,ins);
  }
  void jacFDPre(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi, const double& delta = 1e-6){
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jacFD<0>(J,ins,delta);
  }
  void jacFDPos(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi, const double& delta = 1e-6){
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jacFD<1>(J,ins,delta);
  }
  void jacFDNoi(MXD& J, const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi, const double& delta = 1e-6){
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    this->template _jacFD<2>(J,ins,delta);
  }
  template<int n, int m>
  void setJacBlockPre(MXD& J, const Eigen::Matrix<double,ElementPack<Res...>::template getDim<n>(),ElementPack<Pre...>::template getDim<m>()>& B) const{
    this->template _setJacBlock<0,n,m>(J,B);
  }
  template<int n, int m>
  void setJacBlockPos(MXD& J, const Eigen::Matrix<double,ElementPack<Res...>::template getDim<n>(),ElementPack<Pos...>::template getDim<m>()>& B) const{
    this->template _setJacBlock<1,n,m>(J,B);
  }
  template<int n, int m>
  void setJacBlockNoi(MXD& J, const Eigen::Matrix<double,ElementPack<Res...>::template getDim<n>(),ElementPack<Noi...>::template getDim<m>()>& B) const{
    this->template _setJacBlock<2,n,m>(J,B);
  }
  bool testJacs(const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos, const std::shared_ptr<const StateBase>& noi, const double& delta = 1e-6, const double& th = 1e-6) const{
    const std::array<std::shared_ptr<const StateBase>,3> ins = {pre, pos, noi};
    return this->template _testJacInput<0>(ins,delta,th) &
           this->template _testJacInput<1>(ins,delta,th) &
           this->template _testJacInput<2>(ins,delta,th);
  }

  // Access to definitions
  std::shared_ptr<StateDefinition> resDefinition() const{
    return this->outDefinition_;
  }
  std::shared_ptr<StateDefinition> preDefinition() const{
    return this->inDefinitions_[0];
  }
  std::shared_ptr<StateDefinition> posDefinition() const{
    return this->inDefinitions_[1];
  }
  std::shared_ptr<StateDefinition> noiDefinition() const{
    return this->inDefinitions_[2];
  }

 protected:
  friend mtBase;
  std::shared_ptr<const Meas> meas_;

  // Wrapping from base to user implementation
  void eval(Res&... res, const Pre&... pre, const Pos&... pos, const Noi&... noi) const{
    evalResidual(res...,pre...,pos...,noi...);
  }
  template<int j, typename std::enable_if<(j==0)>::type* = nullptr>
  void jac(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const{
    jacPre(J,pre...,pos...,noi...);
  }
  template<int j, typename std::enable_if<(j==1)>::type* = nullptr>
  void jac(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const{
    jacPos(J,pre...,pos...,noi...);
  }
  template<int j, typename std::enable_if<(j==2)>::type* = nullptr>
  void jac(MXD& J, const Pre&... pre, const Pos&... pos, const Noi&... noi) const{
    jacNoi(J,pre...,pos...,noi...);
  }

};

}

#endif /* GIF_BINARYRESIDUAL_HPP_ */
