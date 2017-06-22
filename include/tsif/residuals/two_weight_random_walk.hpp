#pragma once

#include "tsif/residual.h"

namespace tsif{

class BoolMeasurement : public ElementVector<Element<bool,0>>{
 public:
  BoolMeasurement() : ElementVector<Element<bool,0>>(true){};
  BoolMeasurement(const bool& b) : ElementVector<Element<bool,0>>(b){};
  const bool& get() const {return Get<0>();};
  bool& get() {return Get<0>();};
};

template<typename... Elements>
using TwoWeightRandomWalkBase = Residual<ElementVector<Element<Vec<Elements::kDim>,Elements::kI>...>,
                                ElementVector<Elements...>,
                                ElementVector<Elements...>,
                                BoolMeasurement>;

template<typename... Elements>
class TwoWeightRandomWalk: public TwoWeightRandomWalkBase<Elements...>{
 public:
  typedef ElementVector<Elements...> MyElementVector;
  typedef TwoWeightRandomWalkBase<Elements...> Base;
  using Base::dt_;
  using Base::meas_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  TwoWeightRandomWalk(): Base(true,true,true), w_true_(1.), w_false_(1.) {}
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    _EvalRes(out,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C < MyElementVector::kN)>::type* = nullptr>
  int _EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    const int I = MyElementVector::template GetId<C>();
    cur.template GetElement<I>().Boxminus(pre.template GetElement<I>(),out.template Get<I>());
    _EvalRes<C+1>(out,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C >= MyElementVector::kN)>::type* = nullptr>
  int _EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    _JacPre(J,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C < MyElementVector::kN)>::type* = nullptr>
  int _JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const int I = MyElementVector::template GetId<C>();
    this->template SetJacPre<I,I>(J,pre,cur.template GetElement<I>().BoxminusJacRef(pre.template GetElement<I>()));
    _JacPre<C+1>(J,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C >= MyElementVector::kN)>::type* = nullptr>
  int _JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    _JacCur(J,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C < MyElementVector::kN)>::type* = nullptr>
  int _JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    const int I = MyElementVector::template GetId<C>();
    this->template SetJacCur<I,I>(J,cur,cur.template GetElement<I>().BoxminusJacInp(pre.template GetElement<I>()));
    _JacCur<C+1>(J,pre,cur);
    return 0;
  }
  template<int C=0, typename std::enable_if<(C >= MyElementVector::kN)>::type* = nullptr>
  int _JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  virtual double GetWeight(){
    if (meas_->get()) return w_true_/sqrt(dt_);
    else return w_false_/sqrt(dt_);
  }
  double w_true_;
  double w_false_;
};

} // namespace tsif
