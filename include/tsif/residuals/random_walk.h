#ifndef TSIF_RANDOM_WALK_H_
#define TSIF_RANDOM_WALK_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<typename... Elements>
using RandomWalkBase = Residual<ElementVector<Element<Vec<Elements::kDim>,Elements::kI>...>,
                                ElementVector<Elements...>,
                                ElementVector<Elements...>,
                                MeasEmpty>;

template<typename... Elements>
class RandomWalk: public RandomWalkBase<Elements...>{
 public:
  typedef ElementVector<Elements...> MyElementVector;
  typedef RandomWalkBase<Elements...> Base;
  using Base::dt_;
  using Base::w_;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  RandomWalk(): Base(true,true,true){}
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
    return w_/sqrt(dt_);
  }
};

} // namespace tsif

#endif  // TSIF_RANDOM_WALK_H_
