#ifndef GIF_TRANSFORMATION_HPP_
#define GIF_TRANSFORMATION_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/model.h"

namespace GIF {

template<typename PackOut, typename PackIn>
class Transformation;

template<typename ... Out, typename ... In>
class Transformation<ElementPack<Out...>, ElementPack<In...>> : public Model<
    Transformation<ElementPack<Out...>, ElementPack<In...>>,
    ElementPack<Out...>, ElementPack<In...>> {
 public:
  using mtBase = Model<Transformation<ElementPack<Out...>,ElementPack<In...>>, ElementPack<Out...>,
      ElementPack<In...>>;
  typedef Transformation<ElementPack<Out...>, ElementPack<In...>> mtTransformation;
  Transformation(const std::array<std::string, mtBase::n_>& namesOut,
                 const std::array<std::string, mtBase::m_>& namesIn)
      : mtBase(namesOut, std::forward_as_tuple(namesIn)),
        J_((int) ElementPack<Out...>::d_, (int) ElementPack<In...>::d_) {
  }

  virtual ~Transformation() {
  }


  // User implementations
  virtual void Transform(Out&... outs, const In&... ins) const = 0;
  virtual void JacTransform(MatX& J, const In&... ins) const = 0;

  // Wrapping from user interface to base
  void JacFD(MatX& J, const ElementVectorBase& in, const double delta) {
    const std::array<const ElementVectorBase*, 1> ins = {&in};
    this->template JacFDImpl<0>(J, ins, delta);
  }

  void TransformState(ElementVectorBase* out, const ElementVectorBase& in) {
    const std::array<const ElementVectorBase*, 1> ins = {&in};
    this->template EvalWrapper(out, ins);
  }

  void TransformCovMat(MatX& outputCov, const ElementVectorBase& in, const MatX& inputCov) {
    const std::array<const ElementVectorBase*, 1> ins = {&in};
    this->template JacWrapper<0>(J_, ins);
    outputCov = J_ * inputCov * J_.transpose();
  }

  template<int n, int m>
  void SetJacBlock(
        MatX& J,
        const Eigen::Matrix<double, ElementPack<Out...>::template _GetStateDimension<n>(),
        ElementPack<In...>::template _GetStateDimension<m>()>& B) const {
    this->template SetJacBlockImpl<0, n, m>(J, B);
  }

  bool JacTest(const ElementVectorBase& in, const double delta, const double th) const {
    const std::array<const ElementVectorBase*, 1> ins = {&in};
    return this->template JacTestImpl<0>(ins, delta, th);
  }

  // Access to definitions
  std::shared_ptr<ElementVectorDefinition> OutputDefinition() {
    return this->outDefinition_;
  }
  std::shared_ptr<ElementVectorDefinition> InputDefinition() {
    return this->inDefinitions_[0];
  }

 protected:
  friend mtBase;

  // Wrapping from base to user implementation
  void Eval(Out&... outs, const In&... ins) const {
    Transform(outs..., ins...);
  }

  template<int j, typename std::enable_if<(j == 0)>::type* = nullptr>
  void Jac(MatX& J, const In&... ins) const {
    JacTransform(J, ins...);
  }

  MatX J_;
};

}

#endif /* GIF_TRANSFORMATION_HPP_ */
