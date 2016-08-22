/*
 * Prediction.hpp
 *
 *  Created on: Aug 22, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_PREDICTION_HPP_
#define GIF_PREDICTION_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/BinaryResidual.hpp"

namespace GIF{

template<typename PackSta, typename PackNoi, typename Meas>
class Prediction;

template<typename... Sta, typename... Noi, typename Meas>
class Prediction<ElementPack<Sta...>,ElementPack<Noi...>,Meas>:
    public BinaryResidual<ElementPack<Sta...>,ElementPack<Sta...>,ElementPack<Sta...>,ElementPack<Noi...>,Meas>{
 public:
  typedef Prediction<ElementPack<Sta...>,ElementPack<Noi...>,Meas> mtPrediction;
  typedef BinaryResidual<ElementPack<Sta...>,ElementPack<Sta...>,ElementPack<Sta...>,ElementPack<Noi...>,Meas> mtBinaryRedidual;
  Prediction(const std::array<std::string,ElementPack<Sta...>::n_>& namesSta,
             const std::array<std::string,ElementPack<Noi...>::n_>& namesNoi):
               mtBinaryRedidual(namesSta, namesSta, namesSta, namesNoi, false, true, true){};
  virtual ~Prediction(){};

  // User implementations
  virtual void evalPredictionImpl(Sta&... pos, const Sta&... pre, const Noi&... noi) const = 0;
  virtual void jacPrePredictionImpl(MXD& J, const Sta&... pre, const Noi&... noi) const  = 0;
  virtual void jacNoiPredictionImpl(MXD& J, const Sta&... pre, const Noi&... noi) const  = 0;

 protected:
  // Wrapping from BinaryResidual to Prediction implementation
  void evalResidualImpl(Sta&... res, const Sta&... pre, const Sta&... pos, const Noi&... noi) const{
    // First set innovation to prediction
    evalPredictionImpl(res...,pre...,noi...);
    // Then substract the posterior (inplace)
    computeInnovation(res...,pos...);
  }
  void jacPreImpl(MXD& J, const Sta&... pre, const Sta&... pos, const Noi&... noi) const{
    jacPrePredictionImpl(J,pre...,noi...);
  }
  void jacPosImpl(MXD& J, const Sta&... pre, const Sta&... pos, const Noi&... noi) const{
    J.setIdentity();
    J = -J;
  }
  void jacNoiImpl(MXD& J, const Sta&... pre, const Sta&... pos, const Noi&... noi) const{
    jacNoiPredictionImpl(J,pre...,noi...);
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... res, const Sta&... pos) const{
    typedef typename std::tuple_element<i,typename ElementPack<Sta...>::mtTuple>::type mtElementType;
    Eigen::Matrix<double,ElementTraits<mtElementType>::d_,1> vec;
    ElementTraits<mtElementType>::boxminus(std::get<i>(std::forward_as_tuple(res...)),std::get<i>(std::forward_as_tuple(pos...)),vec);
    ElementTraits<mtElementType>::boxplus(ElementTraits<mtElementType>::identity(),vec,std::get<i>(std::forward_as_tuple(res...))); // TODO: make more efficient
    computeInnovation<i+1>(res...,pos...);
  }
  template<int i = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... res, const Sta&... pos) const{}
};

}

#endif /* GIF_PREDICTION_HPP_ */
