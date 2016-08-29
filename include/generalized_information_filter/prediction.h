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
        prediction_(new ElementVector(this->curDefinition())) {
  }

  virtual ~Prediction() {}

  // User implementations
  virtual void predict(Sta&... cur, const Sta&... pre, const Noi&... noi) const = 0;
  virtual void predictJacPre(MXD& J, const Sta&... pre, const Noi&... noi) const = 0;
  virtual void predictJacNoi(MXD& J, const Sta&... pre, const Noi&... noi) const = 0;

 protected:
  template<typename ... Ts, typename std::enable_if<(sizeof...(Ts)<ElementPack<Sta...>::n_)>::type* = nullptr>
  void _predict(const std::shared_ptr<ElementVectorBase>& cur,
                           const Sta&... pre, const Noi&... noi,
                           Ts&... elements) const {
    assert(cur->MatchesDefinition(this->curDefinition()));
    static constexpr int innerIndex = sizeof...(Ts);
    typedef typename ElementPack<Sta...>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
    _predict(cur, pre..., noi..., elements...,
                        std::dynamic_pointer_cast<Element<mtElementType>>(
                            cur->GetElement(innerIndex))->get());
  }

  template<typename... Ts, typename std::enable_if<(sizeof...(Ts)==ElementPack<Sta...>::n_)>::type* = nullptr>
  void _predict(const std::shared_ptr<ElementVectorBase>& cur,
                           const Sta&... pre, const Noi&... noi,
                           Ts&... elements) const {
    predict(elements..., pre..., noi...);
  }

  // Wrapping from BinaryResidual to Prediction implementation
  void eval(Sta&... inn, const Sta&... pre, const Sta&... cur,
                        const Noi&... noi) const {
    // First compute prediction
    _predict(prediction_, pre..., noi...);
    // Then evaluate difference to posterior
    computeInnovation(inn...,cur...,prediction_);
  }
  void jacPre(MXD& J, const Sta&... pre, const Sta&... cur,
                  const Noi&... noi) const {
    predictJacPre(J,pre...,noi...);
  }
  void jacCur(MXD& J, const Sta&... pre, const Sta&... cur,
                  const Noi&... noi) const {
    J.setZero();
    _predict(prediction_, pre..., noi...);
    computeCurJacobian(J,prediction_,cur...);
  }
  void jacNoi(MXD& J, const Sta&... pre, const Sta&... cur,
                  const Noi&... noi) const {
    predictJacNoi(J,pre...,noi...);
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... inn, const Sta&... cur,
      const std::shared_ptr<const ElementVectorBase>& prediction) const {
    typedef typename std::tuple_element<i,
        typename ElementPack<Sta...>::Tuple>::type mtElementType;

    // inn = I+(pred-cur)
    Eigen::Matrix<double,ElementTraits<mtElementType>::d_,1> vec;
    ElementTraits<mtElementType>::boxminus(
        std::dynamic_pointer_cast<const Element<mtElementType>>(
            prediction->GetElement(i))->get(),std::get<i>(
                std::forward_as_tuple(cur...)),vec);
    // TODO: make more efficient (could be done directly on boxminus, but then
    // jacobian becomes more annoying)
    ElementTraits<mtElementType>::boxplus(
        ElementTraits<mtElementType>::identity(),vec,std::get<i>(
            std::forward_as_tuple(inn...)));
    computeInnovation<i+1>(inn...,cur...,prediction);
  }
  template<int i = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  inline void computeInnovation(Sta&... inn, const Sta&... cur,
      const std::shared_ptr<const ElementVectorBase>& prediction) const {}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void computeCurJacobian(MXD& J, const std::shared_ptr<ElementVectorBase>& prediction,
                          const Sta&... cur) const{
    assert(prediction->MatchesDefinition(this->curDefinition()));
    typedef typename std::tuple_element<i,typename
        ElementPack<Sta...>::Tuple>::type mtElementType;
    Eigen::Matrix<double,ElementTraits<mtElementType>::d_,1> vec;
    ElementTraits<mtElementType>::boxminus(
        std::dynamic_pointer_cast<const Element<mtElementType>>(
            prediction->GetElement(i))->get(),std::get<i>(
                std::forward_as_tuple(cur...)),vec);
    J.template block<ElementTraits<mtElementType>::d_,
                     ElementTraits<mtElementType>::d_>(j,j) =
      ElementTraits<mtElementType>::boxplusJacVec(
          ElementTraits<mtElementType>::identity(),vec) *
      ElementTraits<mtElementType>::boxminusJacRef(
          std::dynamic_pointer_cast<const Element<mtElementType>>(
              prediction->GetElement(i))->get(),std::get<i>(
                  std::forward_as_tuple(cur...)));
    computeCurJacobian<i+1,j+ElementTraits<mtElementType>::d_>(J,prediction,cur...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void computeCurJacobian(MXD& J, const std::shared_ptr<ElementVectorBase>& prediction,
                          const Sta&... cur) const{}

 protected:
  std::shared_ptr<ElementVector> prediction_;
};

}

#endif /* GIF_PREDICTION_HPP_ */
