#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include <list>
#include <initializer_list>

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/element-vector.h"
#include "generalized_information_filter/element-vector-definition.h"

namespace GIF {

/*! \brief Template helper for summing the number of elements.
 */
template<typename ... InPacks>
struct TH_pack_size;

/*! \brief Template helper for getting the index of the first element of a Pack
 */
template<int i, typename ... Packs>
struct TH_pack_index;

/*! \brief Model
 *         A class for implementing the mapping from multiple inputs (ElementVectors) to a single
 *         output (ElementVector). It also provides a finite difference implementation and testing
 *          functionality.
 */
template<typename Derived, typename OutPack, typename ... InPacks>
class Model {
 public:
  template<int InIndex>
  using InPack = typename std::tuple_element<InIndex,std::tuple<InPacks...>>::type;
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  static constexpr int N_ = sizeof...(InPacks);
  typedef Model<Derived,OutPack,InPacks...> mtBase;

  Model(const std::array<std::string,n_>& namesOut,
        const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn);

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<OutPack::n_)>::type* = nullptr>
  inline void EvalWrapper(ElementVectorBase* out,
                          const std::array<const ElementVectorBase*,N_>& ins,
                          Ps&... elements) const;

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps) >= OutPack::n_
               & sizeof...(Ps) <  OutPack::n_ + TH_pack_size<InPacks...>::n_)>::type* = nullptr>
  inline void EvalWrapper(ElementVectorBase* out,
                          const std::array<const ElementVectorBase*,N_>& ins,
                          Ps&... elements) const;

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps) == TH_pack_size<InPacks...>::n_
                                                   +OutPack::n_)>::type* = nullptr>
  inline void EvalWrapper(ElementVectorBase* out,
                          const std::array<const ElementVectorBase*,N_>& ins,
                          Ps&... elements) const;

  template<int InIndex, typename... Ps, typename std::enable_if<
      (sizeof...(Ps)<TH_pack_size<InPacks...>::n_)>::type* = nullptr>
  inline void JacWrapper(MatX& J, const std::array<const ElementVectorBase*,N_>& ins,
                   Ps&... elements) const;

  template<int InIndex, typename... Ps, typename std::enable_if<
      (sizeof...(Ps)==TH_pack_size<InPacks...>::n_)>::type* = nullptr>
  inline void JacWrapper(MatX& J, const std::array<const ElementVectorBase*,N_>& ins,
                   Ps&... elements) const;

  template<int InIndex>
  void JacFDImpl(MatX& J, const std::array<const ElementVectorBase*,N_>& ins,
              const double delta) const;

  template<int InIndex>
  bool JacTestImpl(const std::array<const ElementVectorBase*,N_>& ins,
                     const double delta,
                     const double th) const;

  template<int InIndex>
  bool JacTestImpl(const double delta, const double th) const;

  template<int InIndex, int n, int m>
  MatRef<OutPack::template _GetStateDimension<n>(),
  InPack<InIndex>::template _GetStateDimension<m>()> GetJacBlockImpl(MatRefX J) const;

 protected:
  ElementVectorDefinition::Ptr outDefinition_;
  std::array<ElementVectorDefinition::Ptr,N_> inDefinitions_;

  template<int i = 0, typename std::enable_if<(i<sizeof...(InPacks))>::type* = nullptr>
  void MakeInDefinitions(const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn);

  template<int i = 0, typename std::enable_if<(i==sizeof...(InPacks))>::type* = nullptr>
  void MakeInDefinitions(const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn);
};

// ==================== Implementation ==================== //
template<typename Derived, typename OutPack, typename ... InPacks>
Model<Derived,OutPack,InPacks...>::Model(
      const std::array<std::string,n_>& namesOut,
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
  outDefinition_.reset(new OutPack(namesOut));
  MakeInDefinitions(namesIn);
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<OutPack::n_)>::type*>
inline void Model<Derived,OutPack,InPacks...>::EvalWrapper(
      ElementVectorBase* out,
      const std::array<const ElementVectorBase*,N_>& ins,
      Ps&... elements) const{
  DLOG_IF(FATAL, !out->MatchesDefinition(*outDefinition_)) <<
      "Element vector definition mismatch";
  static constexpr int inner_index = sizeof...(Ps);
  typedef typename OutPack::Tuple Tuple;
  typedef typename std::tuple_element<inner_index,Tuple>::type ElementType;
  EvalWrapper(out, ins, elements..., out->template GetValue<ElementType>(inner_index));
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<typename... Ps,typename std::enable_if<
    (sizeof...(Ps) >= OutPack::n_ & sizeof...(Ps)<OutPack::n_+TH_pack_size<InPacks...>::n_)>::type*>
inline void Model<Derived,OutPack,InPacks...>::EvalWrapper(
      ElementVectorBase* out,
      const std::array<const ElementVectorBase*,N_>& ins,
      Ps&... elements) const{
  static constexpr int outerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::GetOuter();
  static constexpr int innerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::GetInner();
  static_assert(outerIndex < N_, "Indexing Error");
  DLOG_IF(FATAL, !ins.at(outerIndex)->MatchesDefinition(*inDefinitions_[outerIndex])) <<
      "Element vector definition mismatch";
  typedef typename InPack<outerIndex>::Tuple Tuple;
  typedef typename std::tuple_element<innerIndex,Tuple>::type ElementType;
  EvalWrapper(out, ins, elements...,
              ins.at(outerIndex)->template GetValue<ElementType>(innerIndex));
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<typename... Ps, typename std::enable_if<
    (sizeof...(Ps)==TH_pack_size<InPacks...>::n_+OutPack::n_)>::type*>
inline void Model<Derived,OutPack,InPacks...>::EvalWrapper(
      ElementVectorBase* out,
      const std::array<const ElementVectorBase*,N_>& ins,
      Ps&... elements) const{
  static_cast<const Derived&>(*this).Eval(elements...);
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex, typename... Ps, typename std::enable_if<
    (sizeof...(Ps)<TH_pack_size<InPacks...>::n_)>::type*>
inline void Model<Derived,OutPack,InPacks...>::JacWrapper(MatX& J,
                 const std::array<const ElementVectorBase*,N_>& ins,
                 Ps&... elements) const{
  static_assert(InIndex < N_, "Indexing Error");
  static constexpr int outerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::GetOuter();
  static constexpr int innerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::GetInner();
  static_assert(outerIndex < N_, "Indexing Error");
  DLOG_IF(FATAL, !ins.at(outerIndex)->MatchesDefinition(*inDefinitions_[outerIndex])) <<
      "Element vector definition mismatch";
  DLOG_IF(FATAL, J.cols() != inDefinitions_[InIndex]->GetDim()
               | J.rows() != outDefinition_->GetDim()) <<
      "Jacobian Dimension mismatch";
  typedef typename InPack<outerIndex>::Tuple Tuple;
  typedef typename std::tuple_element<innerIndex,Tuple>::type mtElementType;
  JacWrapper<InIndex>(J, ins, elements...,
                ins.at(outerIndex)->template GetValue<mtElementType>(innerIndex));
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex, typename... Ps, typename std::enable_if<
    (sizeof...(Ps)==TH_pack_size<InPacks...>::n_)>::type*>
inline void Model<Derived,OutPack,InPacks...>::JacWrapper(MatX& J,
                 const std::array<const ElementVectorBase*,N_>& ins,
                 Ps&... elements) const{
  static_assert(InIndex<N_,"No such Jacobian!");
  static_cast<const Derived&>(*this).template Jac<InIndex>(J,elements...);
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex>
void Model<Derived,OutPack,InPacks...>::JacFDImpl(MatX& J,
            const std::array<const ElementVectorBase*,N_>& ins,
            const double delta) const{
  static_assert(InIndex < N_, "Indexing Error");
  ElementVector stateDis(inDefinitions_[InIndex]);
  ElementVector outRef(outDefinition_);
  ElementVector outDis(outDefinition_);
  J.resize(outRef.GetDim(),stateDis.GetDim());
  J.setZero();
  stateDis = *ins[InIndex];
  std::array<const ElementVectorBase*,N_> inDis = ins;
  inDis[InIndex] = &stateDis;
  EvalWrapper(&outRef,inDis);
  VecX difIn(stateDis.GetDim());
  VecX difOut(outRef.GetDim());
  for(int i=0; i<stateDis.GetDim(); i++){
    difIn.setZero();
    difIn(i) = delta;
    ins[InIndex]->BoxPlus(difIn,&stateDis);
    EvalWrapper(&outDis,inDis);
    outDis.BoxMinus(outRef,difOut);
    J.col(i) = difOut/delta;
  }
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex>
bool Model<Derived,OutPack,InPacks...>::JacTestImpl(
      const std::array<const ElementVectorBase*,N_>& ins,
      const double delta,
      const double th) const {
  static_assert(InIndex < N_, "Indexing Error");
  if(OutPack::kDim <= 0 || InPack<InIndex>::kDim <= 0){
    return true;
  }
  MatX J((int)OutPack::kDim,(int)InPack<InIndex>::kDim);
  MatX J_FD((int)OutPack::kDim,(int)InPack<InIndex>::kDim);
  JacWrapper<InIndex>(J,ins);
  JacFDImpl<InIndex>(J_FD,ins,delta);
  typename MatX::Index maxRow, maxCol = 0;
  const double r = (J-J_FD).array().abs().maxCoeff(&maxRow, &maxCol);
  if(r>th){
    std::string outName = outDefinition_->GetName(outDefinition_->GetOuter(maxRow));
    std::string inName = inDefinitions_[InIndex]->GetName(ins[InIndex]->GetOuter(maxCol));
    LOG(ERROR) << "==== Model jacInput (" << InIndex << ") Test failed: " << r
               << " is larger than " << th << " at row "
               << maxRow << "("<< outName << "." << outDefinition_->GetInner(maxRow)
               << ") and col " << maxCol << "("<< inName << "."
               << ins[InIndex]->GetInner(maxCol) << ") ====" << std::endl
               << "  " << J(maxRow,maxCol) << "  " << J_FD(maxRow,maxCol)
               << std::endl;
    return false;
  } else {
    LOG(INFO) << "==== Test successful (" << r << ") ====" << std::endl;
    return true;
  }
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex>
bool Model<Derived,OutPack,InPacks...>::JacTestImpl(const double delta, const double th) const{
  std::array<ElementVectorBase::Ptr,N_> ins;
  std::array<const ElementVectorBase*,N_> insRawPtr;
  for(int i=0;i<N_;i++){
    ins[i].reset(new ElementVector(inDefinitions_[i]));
    ins[i]->SetRandom();
    insRawPtr[i] = ins[i].get();
  }
  JacTestImpl<InIndex>(insRawPtr,delta,th);
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int InIndex, int n, int m>
MatRef<OutPack::template _GetStateDimension<n>(),
std::tuple_element<InIndex,std::tuple<InPacks...>>::type::template _GetStateDimension<m>()>
Model<Derived,OutPack,InPacks...>::GetJacBlockImpl(MatRefX J) const{
  return J.block<OutPack::template _GetStateDimension<n>(),
      InPack<InIndex>::template _GetStateDimension<m>()>(
          OutPack::template _GetStartIndex<n>(),
          InPack<InIndex>::template _GetStartIndex<m>());
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int i, typename std::enable_if<(i<sizeof...(InPacks))>::type*>
void Model<Derived,OutPack,InPacks...>::MakeInDefinitions(
    const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
  typedef typename std::tuple_element<i,std::tuple<InPacks...>>::type ElementPack;
  static_assert(i < N_, "Indexing Error");
  inDefinitions_[i].reset(new ElementPack(std::get<i>(namesIn)));
  MakeInDefinitions<i+1>(namesIn);
}

template<typename Derived, typename OutPack, typename ... InPacks>
template<int i, typename std::enable_if<(i==sizeof...(InPacks))>::type*>
void Model<Derived,OutPack,InPacks...>::MakeInDefinitions(
    const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {}

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
