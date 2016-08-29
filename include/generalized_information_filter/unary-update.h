#ifndef GIF_UNARYUPDATE_HPP_
#define GIF_UNARYUPDATE_HPP_

#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/common.h"

namespace GIF {

template<typename PackInn, typename PackSta, typename PackNoi, typename Meas>
class UnaryUpdate;

template<typename ... Inn, typename ... Sta, typename ... Noi, typename Meas>
class UnaryUpdate<ElementPack<Inn...>, ElementPack<Sta...>, ElementPack<Noi...>,
    Meas> : public BinaryResidual<ElementPack<Inn...>, ElementPack<>,
    ElementPack<Sta...>, ElementPack<Noi...>, Meas> {
 public:
  typedef UnaryUpdate<ElementPack<Inn...>, ElementPack<Sta...>,
      ElementPack<Noi...>, Meas> mtUnaryUpdate;
  typedef BinaryResidual<ElementPack<Inn...>, ElementPack<>,
      ElementPack<Sta...>, ElementPack<Noi...>, Meas> mtBinaryRedidual;
  UnaryUpdate(const std::array<std::string, ElementPack<Inn...>::n_>& namesInn,
              const std::array<std::string, ElementPack<Sta...>::n_>& namesSta,
              const std::array<std::string, ElementPack<Noi...>::n_>& namesNoi)
      : mtBinaryRedidual(namesInn, std::array<std::string, 0>(), namesSta,
                         namesNoi, true, false, false) {
  }

  virtual ~UnaryUpdate() {
  }

 protected:
  // Set Jacobian w.r.t. previous state to zero
  void jacPre(MXD& J, const Sta&... sta, const Noi&... noi) const {
    J.setZero();
  }
};

}
#endif /* GIF_UNARYUPDATE_HPP_ */
