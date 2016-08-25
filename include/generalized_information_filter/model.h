#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include <list>
#include <initializer_list>

#include "element-vector-definition.h"
#include "generalized_information_filter/common.h"

namespace GIF {

template<typename ... InPacks>
struct TH_pack_size;

template<int i, typename ... Packs>
struct TH_pack_index;

template<typename Derived, typename OutPack, typename ... InPacks>
class Model {
 public:
  template<int j>
  using InPack = typename std::tuple_element<j,std::tuple<InPacks...>>::type;
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  static constexpr int N_ = sizeof...(InPacks);
  typedef Model<Derived,OutPack,InPacks...> mtBase;

  Model(const std::array<std::string,n_>& namesOut,
        const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    outDefinition_.reset(new OutPack(namesOut));
    MakeInDefinitons(namesIn);
  }

  template<typename... Ps,
           typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  inline void _eval(const SP<ElementVectorBase>& out,
                    const std::array<SP<const ElementVectorBase>,N_>& ins,
                    Ps&... elements) const{
    assert(out->MatchesDefinition(outDefinition_));
    static constexpr int inner_index = sizeof...(Ps);
    typedef typename OutPack::Tuple Tuple;
    typedef typename std::tuple_element<inner_index,Tuple>::type ElementType;
    _eval(out, ins, elements...,
          std::dynamic_pointer_cast<Element<ElementType>>(
              out->GetElement(inner_index))->get());
  }

  template<typename... Ps,
           typename std::enable_if<(sizeof...(Ps)>=n_
                                   & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  inline void _eval(const SP<ElementVectorBase>& out,
                    const std::array<SP<const ElementVectorBase>,N_>& ins,
                    Ps&... elements) const{
    static constexpr int outerIndex =
        TH_pack_index<sizeof...(Ps)-n_,InPacks...>::GetOuter();
    static constexpr int innerIndex =
        TH_pack_index<sizeof...(Ps)-n_,InPacks...>::GetInner();

    assert(ins.at(outerIndex)->MatchesDefinition(inDefinitions_[outerIndex]));
    typedef typename InPack<outerIndex>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type ElementType;
    _eval(out, ins, elements...,
          std::dynamic_pointer_cast<const Element<ElementType>>(
              ins.at(outerIndex)->GetElement(innerIndex))->get());
  }

  template<typename... Ps,
           typename std::enable_if<(sizeof...(Ps)==m_+n_)>::type* = nullptr>
  inline void _eval(const SP<ElementVectorBase>& out,
                    const std::array<SP<const ElementVectorBase>,N_>& ins,
                    Ps&... elements) const{
    static_cast<const Derived&>(*this).eval(elements...);
  }

  template<int j,
           typename... Ps,
           typename std::enable_if<(sizeof...(Ps)<m_)>::type* = nullptr>
  inline void _jac(Mat<>& J,
                   const std::array<SP<const ElementVectorBase>,N_>& ins,
                   Ps&... elements) const{
    static constexpr int outerIndex =
        TH_pack_index<sizeof...(Ps),InPacks...>::GetOuter();
    static constexpr int innerIndex =
        TH_pack_index<sizeof...(Ps),InPacks...>::GetInner();
    assert(ins.at(outerIndex)->MatchesDefinition(inDefinitions_[outerIndex]));
    typedef typename InPack<outerIndex>::Tuple Tuple;
    typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
    _jac<j>(J, ins, elements...,
            std::dynamic_pointer_cast<const Element<mtElementType>>(
                ins.at(outerIndex)->GetElement(innerIndex))->get());
  }

  template<int j,
           typename... Ps,
           typename std::enable_if<(sizeof...(Ps)==m_)>::type* = nullptr>
  inline void _jac(Mat<>& J,
                   const std::array<SP<const ElementVectorBase>,N_>& ins,
                   Ps&... elements) const{
    static_assert(j<N_,"No such Jacobian!");
    static_cast<const Derived&>(*this).template jac<j>(J,elements...);
  }

  template<int j>
  void _jacFD(Mat<>& J,
              const std::array<SP<const ElementVectorBase>,N_>& ins,
              const double& delta = 1e-8) const{
    SP<ElementVector> stateDis(new ElementVector(inDefinitions_[j]));
    SP<ElementVector> outRef(new ElementVector(outDefinition_));
    SP<ElementVector> outDis(new ElementVector(outDefinition_));
    J.resize(outRef->GetDimension(),stateDis->GetDimension());
    J.setZero();
    *stateDis = *ins[j];
    std::array<SP<const ElementVectorBase>,N_> inDis = ins;
    inDis[j] = stateDis;
    _eval(outRef,inDis);
    Vec<> difIn(stateDis->GetDimension());
    Vec<> difOut(outRef->GetDimension());
    for(int i=0; i<stateDis->GetDimension(); i++){
      difIn.setZero();
      difIn(i) = delta;
      ins[j]->BoxPlus(difIn,stateDis);
      _eval(outDis,inDis);
      outDis->BoxMinus(outRef,difOut);
      J.col(i) = difOut/delta;
    }
  }

  template<int j>
  bool _testJacInput(const std::array<SP<const ElementVectorBase>,N_>& ins,
                     const double& delta = 1e-6,
                     const double& th = 1e-6) const {
    Mat<> J((int)OutPack::d_,(int)InPack<j>::d_);
    Mat<> J_FD((int)OutPack::d_,(int)InPack<j>::d_);
    SP<ElementVector> output(new ElementVector(outDefinition_));
    _jac<j>(J,ins);
    _jacFD<j>(J_FD,ins,delta);
    typename Mat<>::Index maxRow, maxCol = 0;
    const double r = (J-J_FD).array().abs().maxCoeff(&maxRow, &maxCol);
    if(r>th){
      std::string outName = outDefinition_->GetName(output->GetOuter(maxRow));
      std::string inName = inDefinitions_[j]->GetName(ins[j]->GetOuter(maxCol));
      std::cout << "==== Model jacInput (" << j << ") Test failed: " << r
                << " is larger than " << th << " at row "
                << maxRow << "("<< outName << "." << output->GetInner(maxRow)
                << ") and col " << maxCol << "("<< inName << "."
                << ins[j]->GetInner(maxCol) << ") ====" << std::endl;
      std::cout << "  " << J(maxRow,maxCol) << "  " << J_FD(maxRow,maxCol)
                << std::endl;
      return false;
    } else {
      std::cout << "==== Test successful (" << r << ") ====" << std::endl;
      return true;
    }
  }

  template<int j>
  bool _testJacInput(int& s, const double& delta = 1e-6,
                     const double& th = 1e-6) const{
    std::array<SP<const ElementVectorBase>,N_> ins;
    for(int i=0;i<N_;i++){
      SP<ElementVectorBase> randomState(new ElementVector(inDefinitions_[i]));
      randomState->SetRandom(s);
      ins[i] = randomState;
    }
    _testJacInput<j>(ins,delta,th);
  }

  template<int j, int n, int m>
  void _setJacBlock(Mat<>& J,
                    const Mat<OutPack::template _GetStateDimension<n>(),
                    std::tuple_element<j,std::tuple<InPacks...>>::type::
                        template _GetStateDimension<m>()>& B) const {
    J.block<OutPack::template _GetStateDimension<n>(),
            InPack<j>::template _GetStateDimension<m>()>(
                OutPack::template _GetStartIndex<n>(),
                InPack<j>::template _GetStartIndex<m>()) = B;
  }

 protected:
  SP<ElementVectorDefinition> outDefinition_;
  std::array<SP<ElementVectorDefinition>,N_> inDefinitions_;

  template<int i = 0, typename std::enable_if<(i<N_)>::type* = nullptr>
  void MakeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    inDefinitions_[i].reset(
        new typename std::tuple_element<i,std::tuple<InPacks...>>::type(
            std::get<i>(namesIn)));
    MakeInDefinitons<i+1>(namesIn);
  }

  template<int i = 0, typename std::enable_if<(i==N_)>::type* = nullptr>
  void MakeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {}
};

// ==================== Implementation ==================== //
template<typename ... Ts, typename ... InPacks>
struct TH_pack_size<ElementPack<Ts...>, InPacks...> {
  static constexpr int n_ = TH_pack_size<InPacks...>::n_ + sizeof...(Ts);
};
template<typename ... Ts>
struct TH_pack_size<ElementPack<Ts...>> {
  static constexpr int n_ = sizeof...(Ts);
};

template<int i, typename ... Ts, typename ... Packs>
struct TH_pack_index<i, ElementPack<Ts...>, Packs...> {
  template<int j = i,
           typename std::enable_if<(j >= sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetOuter() {
    return TH_pack_index<i-sizeof...(Ts),Packs...>::GetOuter()+1;
  }
  template<int j=i,
           typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetOuter() {
    return 0;}
  template<int j=i,
           typename std::enable_if<(j>=sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetInner() {
    return TH_pack_index<i-sizeof...(Ts),Packs...>::GetInner();}
  template<int j=i,
           typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetInner() {
    return i;
  }
};
template<int i, typename ... Ts>
struct TH_pack_index<i, ElementPack<Ts...>> {
  template<typename std::enable_if<(i < sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetOuter() {return 0;};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static constexpr int GetInner() {return i;};
};

}

#endif /* GIF_MODEL_HPP_ */
