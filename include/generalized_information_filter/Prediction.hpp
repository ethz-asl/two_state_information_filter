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

  template<typename... Ts>
  inline void computeInnovation(Ts&... res, const Ts&... pos) const;

  template<typename T>
  inline void computeInnovation(T& r, const T& p) const{
    Eigen::Matrix<double,ElementTraits<T>::d_,1> vec;
    ElementTraits<T>::boxminus(r,p,vec);
    ElementTraits<T>::boxplus(ElementTraits<T>::identity(),vec,r); // TODO: make more efficient
  }

  template<typename T, typename... Ts>
  inline void computeInnovation(T& r, Ts&... res, const T& p, const Ts&... pos) const{
    computeInnovation(r,p);
    computeInnovation(res...,pos...);
  }
};

}

#endif /* GIF_PREDICTION_HPP_ */
