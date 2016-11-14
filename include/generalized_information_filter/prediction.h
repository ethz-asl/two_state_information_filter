#ifndef GIF_PREDICTION_HPP_
#define GIF_PREDICTION_HPP_

#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/common.h"

namespace GIF {

/*! \brief Prediction
 *         Derives from BinaryResidual but enforces identity Jacobian on current state.
 */
template<typename PackSta, typename PackNoi, typename Meas>
class Prediction;

template<typename ... Sta, typename ... Noi, typename Meas>
class Prediction<ElementPack<Sta...>, ElementPack<Noi...>, Meas> :
    public BinaryResidual<ElementPack<Vec<ElementTraits<Sta>::kDim>...>, ElementPack<Sta...>,
                          ElementPack<Sta...>, ElementPack<Noi...>, Meas> {
 public:
  typedef Prediction<ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtPrediction;
  typedef BinaryResidual<ElementPack<Vec<ElementTraits<Sta>::kDim>...>, ElementPack<Sta...>,
                         ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtBinaryRedidual;
  Prediction(const std::string& name,
             const std::array<std::string, ElementPack<Sta...>::n_>& namesSta,
             const std::array<std::string, ElementPack<Noi...>::n_>& namesNoi)
      : mtBinaryRedidual(name, namesSta, namesSta, namesSta, namesNoi, false, true, true),
        prediction_(this->CurDefinition()) {
  }

  virtual ~Prediction() {}

  // User implementations
  virtual void Predict(Sta&... cur, const Sta&... pre, const Noi&... noi) const = 0;
  virtual void PredictJacPre(MatX& J, const Sta&... pre, const Noi&... noi) const = 0;
  virtual void PredictJacNoi(MatX& J, const Sta&... pre, const Noi&... noi) const = 0;

 protected:
  template<typename ... Ts,
           typename std::enable_if<(sizeof...(Ts)<ElementPack<Sta...>::n_)>::type* = nullptr>
  inline void PredictWrapper(ElementVectorBase* cur,
                      const Sta&... pre, const Noi&... noi, Ts&... elements) const {
    DLOG_IF(FATAL,!cur->MatchesDefinition(*this->CurDefinition())) <<
        "Element vector definition mismatch";
    static constexpr int innerIndex = sizeof...(Ts);
    typedef typename ElementPack<Sta...>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
    PredictWrapper(cur, pre..., noi..., elements...,
                   cur->template GetValue<mtElementType>(innerIndex));
  }

  template<typename... Ts,
           typename std::enable_if<(sizeof...(Ts)==ElementPack<Sta...>::n_)>::type* = nullptr>
  inline void PredictWrapper(ElementVectorBase* cur,
                const Sta&... pre, const Noi&... noi, Ts&... elements) const {
    Predict(elements..., pre..., noi...);
  }

  // Wrapping from BinaryResidual to Prediction implementation
  inline void Eval(Vec<ElementTraits<Sta>::kDim>&... inn,
                   const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    // First compute prediction
    PredictWrapper(&prediction_, pre..., noi...);
    // Then evaluate difference to posterior (inn = prediction - cur)
    ComputeInnovation(inn...,cur...,&prediction_);
  }
  inline void JacPre(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    PredictJacPre(J,pre...,noi...);
    PredictWrapper(&prediction_, pre..., noi...); // TODO: cash or avoid
    ComputePreJacobian(J,&prediction_,cur...);
  }
  inline void JacCur(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    J.setZero();
    PredictWrapper(&prediction_, pre..., noi...); // TODO: cash or avoid
    ComputeCurJacobian(J,&prediction_,cur...);
  }
  inline void JacNoi(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    PredictJacNoi(J,pre...,noi...);
    PredictWrapper(&prediction_, pre..., noi...); // TODO: cash or avoid
    ComputeNoiJacobian(J,&prediction_,cur...);
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  inline void ComputeInnovation(Vec<ElementTraits<Sta>::kDim>&... inn, const Sta&... cur,
                                const ElementVectorBase* prediction) const {
    DLOG_IF(FATAL,!prediction->MatchesDefinition(*this->CurDefinition())) <<
        "Element vector definition mismatch";
    typedef typename std::tuple_element<i, typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;
    Trait::Boxminus(prediction->template GetValue<mtElementType>(i),
                    std::get<i>(std::forward_as_tuple(cur...)),
                    std::get<i>(std::forward_as_tuple(inn...)));
    ComputeInnovation<i+1>(inn...,cur...,prediction);
  }
  template<int i = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  inline void ComputeInnovation(Vec<ElementTraits<Sta>::kDim>&... inn, const Sta&... cur,
                                const ElementVectorBase* prediction) const {}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void ComputeCurJacobian(MatX& J, ElementVectorBase* prediction,
                          const Sta&... cur) const{
    DLOG_IF(FATAL,!prediction->MatchesDefinition(*this->CurDefinition())) <<
        "Element vector definition mismatch";
    typedef typename std::tuple_element<i,typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;
    J.template block<Trait::kDim, Trait::kDim>(j,j) =
      Trait::BoxminusJacRef(prediction->template GetValue<mtElementType>(i),
                            std::get<i>(std::forward_as_tuple(cur...)));
    ComputeCurJacobian<i+1,j+Trait::kDim>(J,prediction,cur...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void ComputeCurJacobian(MatX& J, ElementVectorBase* prediction, const Sta&... cur) const{}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void ComputePreJacobian(MatX& J, ElementVectorBase* prediction,
                          const Sta&... cur) const{
    DLOG_IF(FATAL,!prediction->MatchesDefinition(*this->CurDefinition())) <<
        "Element vector definition mismatch";
    typedef typename std::tuple_element<i,typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;
    if(!Trait::kIsVectorSpace){
      J.template block<Trait::kDim, ElementPack<Sta...>::kDim>(j,0) =
        Trait::BoxminusJacInp(prediction->template GetValue<mtElementType>(i),
                              std::get<i>(std::forward_as_tuple(cur...))) *
                              J.template block<Trait::kDim, ElementPack<Sta...>::kDim>(j,0);
    }
    ComputePreJacobian<i+1,j+Trait::kDim>(J,prediction,cur...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void ComputePreJacobian(MatX& J, ElementVectorBase* prediction, const Sta&... cur) const{}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void ComputeNoiJacobian(MatX& J, ElementVectorBase* prediction,
                          const Sta&... cur) const{
    DLOG_IF(FATAL,!prediction->MatchesDefinition(*this->CurDefinition())) <<
        "Element vector definition mismatch";
    typedef typename std::tuple_element<i,typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;
    if(!Trait::kIsVectorSpace){
      J.template block<Trait::kDim, ElementPack<Noi...>::kDim>(j,0) =
        Trait::BoxminusJacInp(prediction->template GetValue<mtElementType>(i),
                              std::get<i>(std::forward_as_tuple(cur...))) *
                              J.template block<Trait::kDim, ElementPack<Sta...>::kDim>(j,0);
    }
    ComputeNoiJacobian<i+1,j+Trait::kDim>(J,prediction,cur...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void ComputeNoiJacobian(MatX& J, ElementVectorBase* prediction, const Sta&... cur) const{}

 protected:
  mutable ElementVector prediction_;
};

}

#endif /* GIF_PREDICTION_HPP_ */
