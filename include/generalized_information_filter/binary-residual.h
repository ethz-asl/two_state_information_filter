#ifndef GIF_BINARYRESIDUAL_HPP_
#define GIF_BINARYRESIDUAL_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/model.h"

namespace GIF {

/*! \brief Base Class for Binary Residual
 *         This class defines the interface for a residual relating two states.
 */
class BinaryResidualBase {
 public:
  BinaryResidualBase(bool isUnary = false, bool isSplitable = false, bool isMergeable = false)
      : isUnary_(isUnary), isSplitable_(isSplitable), isMergeable_(isMergeable){
  }
  virtual ~BinaryResidualBase() {
  }
  virtual void eval(const ElementVectorBase::Ptr& inn,
                            const ElementVectorBase::CPtr& pre,
                            const ElementVectorBase::CPtr& cur,
                            const ElementVectorBase::CPtr& noi) const = 0;
  virtual void jacPre(MatX& J,
                      const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const = 0;
  virtual void jacCur(MatX& J,
                      const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const = 0;
  virtual void jacNoi(MatX& J,
                      const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const = 0;
  virtual ElementVectorDefinition::Ptr innDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr preDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr curDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr noiDefinition() const = 0;
  virtual void splitMeasurements(
      const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
      const ElementVectorBase::CPtr& in,
      ElementVectorBase::CPtr& out1,
      ElementVectorBase::CPtr& out2) const = 0;
  virtual void mergeMeasurements(
      const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
      const ElementVectorBase::CPtr& in1,
      const ElementVectorBase::CPtr& in2,
      ElementVectorBase::CPtr& out) const = 0;
  virtual void setMeas(const ElementVectorBase::CPtr& meas) = 0;
  virtual const MatX& getR() const = 0;
  virtual MatX& getR() = 0;
  virtual bool testJacs(const ElementVectorBase::CPtr& pre,
                        const ElementVectorBase::CPtr& cur,
                        const ElementVectorBase::CPtr& noi,
                        const double& delta = 1e-6, const double& th = 1e-6) const = 0;
  virtual bool testJacs(int& s, const double& delta = 1e-6, const double& th = 1e-6) = 0;

  const bool isUnary_;
  const bool isSplitable_;
  const bool isMergeable_;
};

/*! \brief Binary Residual
 *         A binary residual has three inputs (the two states (previous and current) + an additional
 *         noise term) and returns the innovation. It stores various properties of the residual.
 */
template<typename PackInn, typename PackPre, typename PackCur, typename PackNoi, typename Meas>
class BinaryResidual;

template<typename ... Inn, typename ... Pre, typename ... Cur, typename ... Noi, typename Meas>
class BinaryResidual<ElementPack<Inn...>, ElementPack<Pre...>,ElementPack<Cur...>,
                     ElementPack<Noi...>, Meas>
        : public Model<BinaryResidual<ElementPack<Inn...>, ElementPack<Pre...>, ElementPack<Cur...>,
          ElementPack<Noi...>, Meas>, ElementPack<Inn...>, ElementPack<Pre...>, ElementPack<Cur...>,
          ElementPack<Noi...>>,
          public BinaryResidualBase {
 public:
  typedef BinaryResidual<ElementPack<Inn...>, ElementPack<Pre...>, ElementPack<Cur...>,
      ElementPack<Noi...>, Meas> mtBinaryRedidual;
  using mtBase = Model<mtBinaryRedidual,ElementPack<Inn...>,ElementPack<Pre...>,ElementPack<Cur...>,
      ElementPack<Noi...>>;
  BinaryResidual(const std::array<std::string, ElementPack<Inn...>::n_>& namesInn,
                 const std::array<std::string, ElementPack<Pre...>::n_>& namesPre,
                 const std::array<std::string, ElementPack<Cur...>::n_>& namesCur,
                 const std::array<std::string, ElementPack<Noi...>::n_>& namesNoi,
                 bool isUnary = false, bool isSplitable = false, bool isMergeable = false)
        : mtBase(namesInn, std::forward_as_tuple(namesPre, namesCur, namesNoi)),
          BinaryResidualBase(isUnary, isSplitable, isMergeable),
          meas_(new Meas()) {
    R_.resize(noiDefinition()->GetStateDimension(), noiDefinition()->GetStateDimension());
    R_.setIdentity();
  }
  virtual ~BinaryResidual() {
  }

  // Set measurement
  void setMeas(const ElementVectorBase::CPtr& meas) {
    meas_ = std::dynamic_pointer_cast<const Meas>(meas);
    if (!meas_) {
      assert("ERROR: Passing wrong measurement type" == 0);
      std::cout << "ERROR: Passing wrong measurement type" << std::endl;
    }
  }

  void splitMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                         const ElementVectorBase::CPtr& in,
                         ElementVectorBase::CPtr& out1,
                         ElementVectorBase::CPtr& out2) const {
    // carefull: in/out1/out2 may point to the same elements
    if (isSplitable_) {
      out1 = in;
      out2 = in;
    } else {
      std::cout << "ERROR: splitting of specific residual not supported/implemented!" << std::endl;
    }
  }

  void mergeMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                         const ElementVectorBase::CPtr& in1,
                         const ElementVectorBase::CPtr& in2,
                         ElementVectorBase::CPtr& out) const {
    if (isMergeable_) {
      ElementVectorBase::Ptr newMeas(new Meas());
      VecX diff(in1->GetDimension());
      in1->BoxMinus(in2, diff);
      in2->BoxPlus(toSec(t1 - t0) / toSec(t2 - t0) * diff, newMeas);
      out = newMeas;
    } else {
      std::cout << "ERROR: merging of specific residual not supported!/implemented" << std::endl;
    }
  }

  // User implementations
  virtual void eval(Inn&... inn, const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void jacPre(MatX& J,    const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void jacCur(MatX& J,    const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void jacNoi(MatX& J,    const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;

  // Wrapping from user interface to base
  void eval(const ElementVectorBase::Ptr& inn,
            const ElementVectorBase::CPtr& pre,
            const ElementVectorBase::CPtr& cur,
            const ElementVectorBase::CPtr& noi) const {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->_eval(inn, ins);
  }

  void jacPre(MatX& J, const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jac<0>(J, ins);
  }
  void jacCur(MatX& J, const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jac<1>(J, ins);
  }

  void jacNoi(MatX& J, const ElementVectorBase::CPtr& pre,
                      const ElementVectorBase::CPtr& cur,
                      const ElementVectorBase::CPtr& noi) const {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jac<2>(J, ins);
  }
  void jacFDPre(MatX& J, const ElementVectorBase::CPtr& pre,
                        const ElementVectorBase::CPtr& cur,
                        const ElementVectorBase::CPtr& noi,
                        const double& delta = 1e-6) {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jacFD<0>(J, ins, delta);
  }
  void jacFDCur(MatX& J, const ElementVectorBase::CPtr& pre,
                        const ElementVectorBase::CPtr& cur,
                        const ElementVectorBase::CPtr& noi,
                        const double& delta = 1e-6) {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jacFD<1>(J, ins, delta);
  }

  void jacFDNoi(MatX& J, const ElementVectorBase::CPtr& pre,
                        const ElementVectorBase::CPtr& cur,
                        const ElementVectorBase::CPtr& noi,
                        const double& delta = 1e-6) {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    this->template _jacFD<2>(J, ins, delta);
  }

  template<int n, int m>
  void setJacBlockPre(
      MatX& J,
      const Mat<ElementPack<Inn...>::template _GetStateDimension<n>(),
                ElementPack<Pre...>::template _GetStateDimension<m>()>& B) const {
    this->template _setJacBlock<0, n, m>(J, B);
  }

  template<int n, int m>
  void setJacBlockCur(
      MatX& J,
      const Mat<ElementPack<Inn...>::template _GetStateDimension<n>(),
                ElementPack<Cur...>::template _GetStateDimension<m>()>& B) const {
    this->template _setJacBlock<1, n, m>(J, B);
  }

  template<int n, int m>
  void setJacBlockNoi(
      MatX& J,
      const Mat<ElementPack<Inn...>::template _GetStateDimension<n>(),
                ElementPack<Noi...>::template _GetStateDimension<m>()>& B) const {
    this->template _setJacBlock<2, n, m>(J, B);
  }

  bool testJacs(const ElementVectorBase::CPtr& pre,
                const ElementVectorBase::CPtr& cur,
                const ElementVectorBase::CPtr& noi,
                const double& delta = 1e-6, const double& th = 1e-6) const {
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    return this->template _testJacInput<0>(ins, delta, th)
         & this->template _testJacInput<1>(ins, delta, th)
         & this->template _testJacInput<2>(ins, delta, th);
  }

  bool testJacs(int& s, const double& delta = 1e-6, const double& th = 1e-6) {
    std::shared_ptr<Meas> meas(new Meas());
    meas->SetRandom(s);
    meas_ = meas;
    ElementVectorBase::Ptr pre(new ElementVector(preDefinition()));
    pre->SetRandom(s);
    ElementVectorBase::Ptr cur(new ElementVector(curDefinition()));
    cur->SetRandom(s);
    ElementVectorBase::Ptr noi(new ElementVector(noiDefinition()));
    noi->SetIdentity();
    const std::array<ElementVectorBase::CPtr, 3> ins = {pre, cur, noi};
    return this->template _testJacInput<0>(ins, delta, th)
         & this->template _testJacInput<1>(ins, delta, th)
         & this->template _testJacInput<2>(ins, delta, th);
  }

  // Access to definitions
  ElementVectorDefinition::Ptr innDefinition() const {
    return this->outDefinition_;
  }
  ElementVectorDefinition::Ptr preDefinition() const {
    return this->inDefinitions_[0];
  }
  ElementVectorDefinition::Ptr curDefinition() const {
    return this->inDefinitions_[1];
  }
  ElementVectorDefinition::Ptr noiDefinition() const {
    return this->inDefinitions_[2];
  }

  // Get Noise Matrix
  const MatX& getR() const {
    return R_;
  }
  MatX& getR() {
    return R_;
  }

 protected:
  friend mtBase;
  std::shared_ptr<const Meas> meas_;
  MatX R_;

  // Wrapping from base (model) to user implementation
  template<int j, typename std::enable_if<(j == 0)>::type* = nullptr>
  void jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    jacPre(J, pre..., cur..., noi...);
  }
  template<int j, typename std::enable_if<(j == 1)>::type* = nullptr>
  void jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    jacCur(J, pre..., cur..., noi...);
  }
  template<int j, typename std::enable_if<(j == 2)>::type* = nullptr>
  void jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    jacNoi(J, pre..., cur..., noi...);
  }

};

}
#endif /* GIF_BINARYRESIDUAL_HPP_ */
