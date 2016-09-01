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
    public BinaryResidual<ElementPack<Sta...>, ElementPack<Sta...>,
                          ElementPack<Sta...>, ElementPack<Noi...>, Meas> {
 public:
  typedef Prediction<ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtPrediction;
  typedef BinaryResidual<ElementPack<Sta...>, ElementPack<Sta...>,
                         ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtBinaryRedidual;
  Prediction(const std::array<std::string, ElementPack<Sta...>::n_>& namesSta,
             const std::array<std::string, ElementPack<Noi...>::n_>& namesNoi)
      : mtBinaryRedidual(namesSta, namesSta, namesSta, namesNoi, false, true, true),
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
  void PredictWrapper(ElementVectorBase* cur,
                      const Sta&... pre, const Noi&... noi, Ts&... elements) const {
    assert(cur->MatchesDefinition(*this->CurDefinition()));
    static constexpr int innerIndex = sizeof...(Ts);
    typedef typename ElementPack<Sta...>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
    PredictWrapper(cur, pre..., noi..., elements...,
                   cur->template GetValue<mtElementType>(innerIndex));
  }

  template<typename... Ts,
           typename std::enable_if<(sizeof...(Ts)==ElementPack<Sta...>::n_)>::type* = nullptr>
  void PredictWrapper(ElementVectorBase* cur,
                const Sta&... pre, const Noi&... noi, Ts&... elements) const {
    Predict(elements..., pre..., noi...);
  }

  // Wrapping from BinaryResidual to Prediction implementation
  void Eval(Sta&... inn, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    // First compute prediction
    PredictWrapper(&prediction_, pre..., noi...);
    // Then evaluate difference to posterior
    ComputeInnovation(inn...,cur...,&prediction_);
  }
  void JacPre(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    PredictJacPre(J,pre...,noi...);
  }
  void JacCur(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    J.setZero();
    PredictWrapper(&prediction_, pre..., noi...); // TODO: cash or avoid
    ComputeCurJacobian(J,&prediction_,cur...);
  }
  void JacNoi(MatX& J, const Sta&... pre, const Sta&... cur, const Noi&... noi) const {
    PredictJacNoi(J,pre...,noi...);
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  inline void ComputeInnovation(Sta&... inn, const Sta&... cur,
                                const ElementVectorBase* prediction) const {
    typedef typename std::tuple_element<i, typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;

    // inn = I+(pred-cur)
    // TODO: make more efficient (could be done directly on Boxminus, but then jacobian becomes more annoying)
    Eigen::Matrix<double,Trait::kDim,1> vec;
    Trait::Boxminus(prediction->template GetValue<mtElementType>(i),
                    std::get<i>(std::forward_as_tuple(cur...)),
                    vec);
    Trait::Boxplus(Trait::Identity(), vec, std::get<i>(std::forward_as_tuple(inn...)));
    ComputeInnovation<i+1>(inn...,cur...,prediction);
  }
  template<int i = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  inline void ComputeInnovation(Sta&... inn, const Sta&... cur,
                                const ElementVectorBase* prediction) const {}

  template<int i = 0, int j = 0, typename std::enable_if<(i<sizeof...(Sta))>::type* = nullptr>
  void ComputeCurJacobian(MatX& J, ElementVectorBase* prediction,
                          const Sta&... cur) const{
    assert(prediction->MatchesDefinition(*this->CurDefinition()));
    typedef typename std::tuple_element<i,typename ElementPack<Sta...>::Tuple>::type mtElementType;
    typedef ElementTraits<mtElementType> Trait;
    Eigen::Matrix<double,Trait::kDim,1> vec;
    Trait::Boxminus(prediction->template GetValue<mtElementType>(i),
                    std::get<i>(std::forward_as_tuple(cur...)),
                    vec);
    J.template block<Trait::kDim, Trait::kDim>(j,j) = Trait::BoxplusJacVec(Trait::Identity(),vec) *
      Trait::BoxminusJacRef(prediction->template GetValue<mtElementType>(i),
                            std::get<i>(std::forward_as_tuple(cur...)));
    ComputeCurJacobian<i+1,j+Trait::kDim>(J,prediction,cur...);
  }
  template<int i = 0, int j = 0, typename std::enable_if<(i>=sizeof...(Sta))>::type* = nullptr>
  void ComputeCurJacobian(MatX& J, ElementVectorBase* prediction, const Sta&... cur) const{}

 protected:
  mutable ElementVector prediction_;
};

}

#endif /* GIF_PREDICTION_HPP_ */
