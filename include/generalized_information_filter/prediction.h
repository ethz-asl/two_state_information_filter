#ifndef GIF_PREDICTION_HPP_
#define GIF_PREDICTION_HPP_

#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/common.h"

namespace GIF {

template<typename PackSta, typename PackNoi, typename Meas>
class Prediction;

template<typename ... Sta, typename ... Noi, typename Meas>
class Prediction<ElementPack<Sta...>, ElementPack<Noi...>, Meas> :
    public BinaryResidual<ElementPack<Sta...>, ElementPack<Sta...>,
        ElementPack<Sta...>, ElementPack<Noi...>, Meas> {
 public:
  typedef Prediction<ElementPack<Sta...>, ElementPack<Noi...>, Meas>
    mtPrediction;
  typedef BinaryResidual<ElementPack<Sta...>, ElementPack<Sta...>,
      ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtBinaryRedidual;
  Prediction(const std::array<std::string, ElementPack<Sta...>::n_>& namesSta,
             const std::array<std::string, ElementPack<Noi...>::n_>& namesNoi)
      : mtBinaryRedidual(namesSta, namesSta, namesSta, namesNoi, false, true,
                         true),
        prediction_(new ElementVector(this->posDefinition())) {
  }

  virtual ~Prediction() {}

  // User implementations
  virtual void evalPredictionImpl(Sta&... pos, const Sta&... pre,
                                  const Noi&... noi) const = 0;
  virtual void jacPrePredictionImpl(MXD& J, const Sta&... pre,
                                    const Noi&... noi) const = 0;
  virtual void jacNoiPredictionImpl(MXD& J, const Sta&... pre,
                                    const Noi&... noi) const = 0;

 protected:
  template<typename ... Ts, typename std::enable_if<(sizeof...(Ts)<ElementPack<Sta...>::n_)>::type* = nullptr>
  void _evalPredictionImpl(const std::shared_ptr<ElementVectorBase>& pos,
                           const Sta&... pre, const Noi&... noi,
                           Ts&... elements) const {
    assert(pos->MatchesDefinition(this->posDefinition()));
    static constexpr int innerIndex = sizeof...(Ts);
    typedef typename ElementPack<Sta...>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
    _evalPredictionImpl(pos, pre..., noi..., elements...,
                        std::dynamic_pointer_cast<Element<mtElementType>>(
                            pos->GetElement(innerIndex))->get());
  }

  template<typename... Ts, typename std::enable_if<(sizeof...(Ts)==ElementPack<Sta...>::n_)>::type* = nullptr>
  void _evalPredictionImpl(const std::shared_ptr<ElementVectorBase>& pos,
                           const Sta&... pre, const Noi&... noi,
                           Ts&... elements) const {
    evalPredictionImpl(elements..., pre..., noi...);
  }

  // Wrapping from BinaryResidual to Prediction implementation
  void evalResidualImpl(Sta&... res, const Sta&... pre, const Sta&... pos,
                        const Noi&... noi) const {
    // First compute prediction
    _evalPredictionImpl(prediction_, pre..., noi...);
    // Then evaluate difference to posterior
    computeInnovation(res...,pos...,prediction_);
  }
  void jacPreImpl(MXD& J, const Sta&... pre, const Sta&... pos,
                  const Noi&... noi) const {
    jacPrePredictionImpl(J,pre...,noi...);
  }
  void jacPosImpl(MXD& J, const Sta&... pre, const Sta&... pos,
                  const Noi&... noi) const {
    J.setZero();
    _evalPredictionImpl(prediction_, pre..., noi...);
    computePosJacobian(J,prediction_,pos...);
  }
  void jacNoiImpl(MXD& J, const Sta&... pre, const Sta&... pos,
                  const Noi&... noi) const {
    jacNoiPredictionImpl(J,pre...,noi...);
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... res, const Sta&... pos,
      const std::shared_ptr<const ElementVectorBase>& prediction) const {
    typedef typename std::tuple_element<i,
        typename ElementPack<Sta...>::Tuple>::type mtElementType;

    // res = I+(pred-pos)
    Eigen::Matrix<double,ElementTraits<mtElementType>::d_,1> vec;
    ElementTraits<mtElementType>::boxminus(
        std::dynamic_pointer_cast<const Element<mtElementType>>(
            prediction->GetElement(i))->get(),std::get<i>(
                std::forward_as_tuple(pos...)),vec);
    // TODO: make more efficient (could be done directly on boxminus, but then
    // jacobian becomes more annoying)
    ElementTraits<mtElementType>::boxplus(
        ElementTraits<mtElementType>::identity(),vec,std::get<i>(
            std::forward_as_tuple(res...)));
    computeInnovation<i+1>(res...,pos...,prediction);
  }
  template<int i = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... res, const Sta&... pos,
      const std::shared_ptr<const ElementVectorBase>& prediction) const {}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void computePosJacobian(MXD& J, const std::shared_ptr<ElementVectorBase>& prediction,
                          const Sta&... pos) const{
    assert(prediction->MatchesDefinition(this->posDefinition()));
    typedef typename std::tuple_element<i,typename
        ElementPack<Sta...>::Tuple>::type mtElementType;
    Eigen::Matrix<double,ElementTraits<mtElementType>::d_,1> vec;
    ElementTraits<mtElementType>::boxminus(
        std::dynamic_pointer_cast<const Element<mtElementType>>(
            prediction->GetElement(i))->get(),std::get<i>(
                std::forward_as_tuple(pos...)),vec);
    J.template block<ElementTraits<mtElementType>::d_,
                     ElementTraits<mtElementType>::d_>(j,j) =
      ElementTraits<mtElementType>::boxplusJacVec(
          ElementTraits<mtElementType>::identity(),vec) *
      ElementTraits<mtElementType>::boxminusJacRef(
          std::dynamic_pointer_cast<const Element<mtElementType>>(
              prediction->GetElement(i))->get(),std::get<i>(
                  std::forward_as_tuple(pos...)));
    computePosJacobian<i+1,j+ElementTraits<mtElementType>::d_>(J,prediction,pos...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void computePosJacobian(MXD& J, const std::shared_ptr<ElementVectorBase>& prediction,
                          const Sta&... pos) const{}

 protected:
  std::shared_ptr<ElementVector> prediction_;
};

}

#endif /* GIF_PREDICTION_HPP_ */
