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
  typedef std::shared_ptr<BinaryResidualBase> Ptr;
  typedef std::shared_ptr<const BinaryResidualBase> CPtr;
  BinaryResidualBase(bool isUnary = false, bool isSplitable = false, bool isMergeable = false)
      : isUnary_(isUnary), isSplitable_(isSplitable), isMergeable_(isMergeable){
    dt_ = 0.1;
  }
  virtual ~BinaryResidualBase() {
  }
  virtual void Eval(ElementVectorBase* inn,
                    const ElementVectorBase& pre,
                    const ElementVectorBase& cur,
                    const ElementVectorBase& noi) const = 0;
  virtual void JacPre(MatX& J,
                      const ElementVectorBase& pre,
                      const ElementVectorBase& cur,
                      const ElementVectorBase& noi) const = 0;
  virtual void JacCur(MatX& J,
                      const ElementVectorBase& pre,
                      const ElementVectorBase& cur,
                      const ElementVectorBase& noi) const = 0;
  virtual void JacNoi(MatX& J,
                      const ElementVectorBase& pre,
                      const ElementVectorBase& cur,
                      const ElementVectorBase& noi) const = 0;
  virtual ElementVectorDefinition::Ptr InnDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr PreDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr CurDefinition() const = 0;
  virtual ElementVectorDefinition::Ptr NoiDefinition() const = 0;
  virtual void SplitMeasurements(
      const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
      const ElementVectorBase::CPtr& in,
      ElementVectorBase::CPtr& out1,
      ElementVectorBase::CPtr& out2) const = 0;
  virtual void MergeMeasurements(
      const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
      const ElementVectorBase::CPtr& in1,
      const ElementVectorBase::CPtr& in2,
      ElementVectorBase::CPtr& out) const = 0;
  virtual void SetMeas(const ElementVectorBase::CPtr& meas) = 0;
  void SetDt(double dt){
    dt_ = dt;
  }
  virtual bool CheckMeasType(const ElementVectorBase::CPtr& meas) const = 0;
  virtual const MatX& GetNoiseCovariance() const = 0;
  virtual MatX& GetNoiseCovariance() = 0;
  virtual bool TestJacs(const ElementVectorBase& pre,
                        const ElementVectorBase& cur,
                        const ElementVectorBase& noi,
                        const double delta, const double th) const = 0;
  virtual bool TestJacs(const double delta, const double th) = 0;

  const bool isUnary_;
  const bool isSplitable_;
  const bool isMergeable_;
  double dt_;
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
                 bool isUnary, bool isSplitable, bool isMergeable)
        : mtBase(namesInn, std::forward_as_tuple(namesPre, namesCur, namesNoi)),
          BinaryResidualBase(isUnary, isSplitable, isMergeable),
          meas_(new Meas()) {
    R_.resize(NoiDefinition()->GetDim(), NoiDefinition()->GetDim());
    R_.setIdentity();
  }
  virtual ~BinaryResidual() {
  }

  // Set measurement
  void SetMeas(const ElementVectorBase::CPtr& meas) {
    meas_ = std::dynamic_pointer_cast<const Meas>(meas);
    DLOG_IF(ERROR, !meas_) << "Passing wrong measurement type";
  }

  bool CheckMeasType(const ElementVectorBase::CPtr& meas) const{
    std::shared_ptr<const Meas> cast_meas = std::dynamic_pointer_cast<const Meas>(meas);
    return cast_meas.get() != nullptr;
  }

  void SplitMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                         const ElementVectorBase::CPtr& in,
                         ElementVectorBase::CPtr& out1,
                         ElementVectorBase::CPtr& out2) const {
    // carefull: in/out1/out2 may point to the same elements
    if (isSplitable_) {
      out1 = in;
      out2 = in;
    } else {
      LOG(ERROR) << "Splitting of specific residual not supported/implemented!";
    }
  }

  void MergeMeasurements(const TimePoint& t0, const TimePoint& t1, const TimePoint& t2,
                         const ElementVectorBase::CPtr& in1,
                         const ElementVectorBase::CPtr& in2,
                         ElementVectorBase::CPtr& out) const {
    if (isMergeable_) {
      ElementVectorBase::Ptr newMeas(new Meas());
      VecX diff(in1->GetDim());
      in1->BoxMinus(*in2, diff);
      in2->BoxPlus(toSec(t1 - t0) / toSec(t2 - t0) * diff, newMeas.get());
      out = newMeas;
    } else {
      LOG(ERROR) << "Merging of specific residual not supported/implemented!";
    }
  }

  // User implementations
  virtual void Eval(Inn&... inn, const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void JacPre(MatX& J,   const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void JacCur(MatX& J,   const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;
  virtual void JacNoi(MatX& J,   const Pre&... pre, const Cur&... cur, const Noi&... noi) const = 0;

  // Interface to base
  inline void Eval(ElementVectorBase* inn,
            const ElementVectorBase& pre,
            const ElementVectorBase& cur,
            const ElementVectorBase& noi) const {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->EvalWrapper(inn, ins);
  }

  inline void JacPre(MatX& J, const ElementVectorBase& pre,
                       const ElementVectorBase& cur,
                       const ElementVectorBase& noi) const {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacWrapper<0>(J, ins);
  }
  inline void JacCur(MatX& J, const ElementVectorBase& pre,
                       const ElementVectorBase& cur,
                       const ElementVectorBase& noi) const {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacWrapper<1>(J, ins);
  }

  inline void JacNoi(MatX& J, const ElementVectorBase& pre,
                       const ElementVectorBase& cur,
                       const ElementVectorBase& noi) const {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacWrapper<2>(J, ins);
  }
  inline void JacFDPre(MatX& J, const ElementVectorBase& pre,
                         const ElementVectorBase& cur,
                         const ElementVectorBase& noi,
                         const double delta) {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacFDImpl<0>(J, ins, delta);
  }
  inline void JacFDCur(MatX& J, const ElementVectorBase& pre,
                         const ElementVectorBase& cur,
                         const ElementVectorBase& noi,
                         const double delta) {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacFDImpl<1>(J, ins, delta);
  }

  inline void JacFDNoi(MatX& J, const ElementVectorBase& pre,
                         const ElementVectorBase& cur,
                         const ElementVectorBase& noi,
                         const double delta) {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    this->template JacFDImpl<2>(J, ins, delta);
  }

  template<int n, int m>
  inline MatRef<ElementPack<Inn...>::template _GetStateDimension<n>(),
         ElementPack<Pre...>::template _GetStateDimension<m>()> GetJacBlockPre(MatRefX J) const {
    return this->template GetJacBlockImpl<0, n, m>(J);
  }

  template<int n, int m>
  inline MatRef<ElementPack<Inn...>::template _GetStateDimension<n>(),
         ElementPack<Cur...>::template _GetStateDimension<m>()> GetJacBlockCur(MatRefX J) const {
    return this->template GetJacBlockImpl<1, n, m>(J);
  }

  template<int n, int m>
  inline MatRef<ElementPack<Inn...>::template _GetStateDimension<n>(),
         ElementPack<Noi...>::template _GetStateDimension<m>()> GetJacBlockNoi(MatRefX J) const {
    return this->template GetJacBlockImpl<2, n, m>(J);
  }

  bool TestJacs(const ElementVectorBase& pre,
                const ElementVectorBase& cur,
                const ElementVectorBase& noi,
                const double delta, const double th) const {
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    return this->template JacTestImpl<0>(ins, delta, th)
         & this->template JacTestImpl<1>(ins, delta, th)
         & this->template JacTestImpl<2>(ins, delta, th);
  }

  bool TestJacs(const double delta, const double th) {
    std::shared_ptr<Meas> meas(new Meas());
    meas->SetRandom();
    meas_ = meas;
    ElementVector pre(PreDefinition());
    pre.SetRandom();
    ElementVector cur(CurDefinition());
    cur.SetRandom();
    ElementVector noi(NoiDefinition());
    noi.SetIdentity();
    const std::array<const ElementVectorBase*, 3> ins = {&pre, &cur, &noi};
    return this->template JacTestImpl<0>(ins, delta, th)
         & this->template JacTestImpl<1>(ins, delta, th)
         & this->template JacTestImpl<2>(ins, delta, th);
  }

  // Access to definitions
  inline ElementVectorDefinition::Ptr InnDefinition() const {
    return this->outDefinition_;
  }
  inline ElementVectorDefinition::Ptr PreDefinition() const {
    return this->inDefinitions_[0];
  }
  inline ElementVectorDefinition::Ptr CurDefinition() const {
    return this->inDefinitions_[1];
  }
  inline ElementVectorDefinition::Ptr NoiDefinition() const {
    return this->inDefinitions_[2];
  }

  // Get Noise Matrix
  inline const MatX& GetNoiseCovariance() const {
    return R_;
  }
  inline MatX& GetNoiseCovariance() {
    return R_;
  }

 protected:
  friend mtBase;
  std::shared_ptr<const Meas> meas_;
  MatX R_;

  // Wrapping from base (model) to user implementation
  template<int j, typename std::enable_if<(j == 0)>::type* = nullptr>
  inline void Jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    JacPre(J, pre..., cur..., noi...);
  }
  template<int j, typename std::enable_if<(j == 1)>::type* = nullptr>
  inline void Jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    JacCur(J, pre..., cur..., noi...);
  }
  template<int j, typename std::enable_if<(j == 2)>::type* = nullptr>
  inline void Jac(MatX& J, const Pre&... pre, const Cur&... cur, const Noi&... noi) const {
    JacNoi(J, pre..., cur..., noi...);
  }

};

}
#endif /* GIF_BINARYRESIDUAL_HPP_ */
