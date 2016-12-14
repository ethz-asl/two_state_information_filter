#ifndef GIF_RWPREDICTION_HPP_
#define GIF_RWPREDICTION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/prediction.h"
#include "generalized_information_filter/measurement.h"

namespace GIF {


/*! \brief Random Walk Prediction.
 *         Predicts the current state based on the previous state assuming a random walk model.
 */
template<typename PackSta>
class RandomWalkPrediction;

template<typename ... Sta>
class RandomWalkPrediction<ElementPack<Sta...>>: public Prediction<ElementPack<Sta...>,
                                        ElementPack<Vec<ElementTraits<Sta>::kDim>...>,
                                        EmptyMeas> {
 public:
  static constexpr int n_ = ElementPack<Sta...>::n_;
  typedef std::array<std::string,n_> StringArr;
  using mtPrediction = Prediction<ElementPack<Sta...>,
                                  ElementPack<Vec<ElementTraits<Sta>::kDim>...>,
                                  EmptyMeas>;
  using mtPrediction::dt_;
  RandomWalkPrediction(const std::string& name, const StringArr& staName, const StringArr& noiName)
      : mtPrediction(name, staName, noiName){
  }
  virtual ~RandomWalkPrediction() {
  }
  void Predict(Sta&... cur, const Sta&... pre, const Vec<ElementTraits<Sta>::kDim>&... noi) const {
    PredictRecursion(std::forward_as_tuple(cur...),
                     std::forward_as_tuple(pre...),
                     std::forward_as_tuple(noi...));
  }
  template<int i=0, typename std::enable_if<i<n_>::type* = nullptr>
  void PredictRecursion(std::tuple<Sta&...> cur,
                        std::tuple<const Sta&...> pre,
                        std::tuple<const Vec<ElementTraits<Sta>::kDim>&...> noi) const {
    ElementTraits<typename std::tuple_element<i,std::tuple<Sta...>>::type>::Boxplus(
        std::get<i>(pre), std::get<i>(noi) * sqrt(dt_), std::get<i>(cur));
    PredictRecursion<i+1>(cur,pre,noi);
  }
  template<int i=0, typename std::enable_if<i>=n_>::type* = nullptr>
  void PredictRecursion(std::tuple<Sta&...> cur,
                        std::tuple<const Sta&...> pre,
                        std::tuple<const Vec<ElementTraits<Sta>::kDim>&...> noi) const {}
  void PredictJacPre(MatX& J,
                     const Sta&... pre, const Vec<ElementTraits<Sta>::kDim>&... noi) const {
    J.setIdentity();
  }
  void PredictJacNoi(MatX& J,
                     const Sta&... pre, const Vec<ElementTraits<Sta>::kDim>&... noi) const {
    J.setIdentity();
    J = sqrt(dt_)*J; // TODO: make more efficient
  }
};

}
#endif /* GIF_RWPREDICTION_HPP_ */
