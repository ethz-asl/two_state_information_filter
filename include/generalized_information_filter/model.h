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
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  static constexpr int N_ = sizeof...(InPacks);
  typedef Model<Derived,OutPack,InPacks...> mtBase;

  Model(const std::array<std::string,n_>& namesOut,
        const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    outDefinition_.reset(new OutPack(namesOut));
    makeInDefinitons(namesIn);
  }

  template<int i = 0, typename std::enable_if<(i<N_)>::type* = nullptr>
  void makeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    inDefinitions_[i].reset(
        new typename std::tuple_element<i,std::tuple<InPacks...>>::type(
            std::get<i>(namesIn)));
    makeInDefinitons<i+1>(namesIn);
  }

  template<int i = 0, typename std::enable_if<(i==N_)>::type* = nullptr>
  void makeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {}

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::shared_ptr<ElementVectorBase>& out,
             const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
             Ps&... elements) const{
    assert(out->matchesDef(outDefinition_));
    static constexpr int innerIndex = sizeof...(Ps);
    typedef typename OutPack::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _eval(out, ins, elements...,
          std::dynamic_pointer_cast<Element<mtElementType>>(
              out->getElement(innerIndex))->get());
  }

  template<typename... Ps,
           typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _eval(const std::shared_ptr<ElementVectorBase>& out,
             const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
             Ps&... elements) const{
    static constexpr int outerIndex =
        TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getOuter();
    static constexpr int innerIndex =
        TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getInner();

    assert(ins.at(outerIndex)->matchesDef(inDefinitions_[outerIndex]));
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _eval(out, ins, elements...,
          std::dynamic_pointer_cast<const Element<mtElementType>>(
              ins.at(outerIndex)->getElement(innerIndex))->get());
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==m_+n_)>::type* = nullptr>
  void _eval(const std::shared_ptr<ElementVectorBase>& out,
             const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
             Ps&... elements) const{
    static_cast<const Derived&>(*this).eval(elements...);
  }

  template<int j, typename... Ps, typename std::enable_if<(sizeof...(Ps)<m_)>::type* = nullptr>
  void _jac(MXD& J, const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
            Ps&... elements) const{
    static constexpr int outerIndex =
        TH_pack_index<sizeof...(Ps),InPacks...>::getOuter();
    static constexpr int innerIndex =
        TH_pack_index<sizeof...(Ps),InPacks...>::getInner();
    assert(ins.at(outerIndex)->matchesDef(inDefinitions_[outerIndex]));
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _jac<j>(J, ins, elements...,
            std::dynamic_pointer_cast<const Element<mtElementType>>(
                ins.at(outerIndex)->getElement(innerIndex))->get());
  }

  template<int j, typename... Ps, typename std::enable_if<(sizeof...(Ps)==m_)>::type* = nullptr>
  void _jac(MXD& J, const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
            Ps&... elements) const{
    static_assert(j<N_,"No such Jacobian!");
    static_cast<const Derived&>(*this).template jac<j>(J,elements...);
  }

  template<int j>
  void _jacFD(MXD& J,
              const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
              const double& delta = 1e-8) const{
    std::shared_ptr<ElementVector> stateDis(new ElementVector(inDefinitions_[j]));
    std::shared_ptr<ElementVector> outRef(new ElementVector(outDefinition_));
    std::shared_ptr<ElementVector> outDis(new ElementVector(outDefinition_));
    J.resize(outRef->getDim(),stateDis->getDim());
    J.setZero();
    *stateDis = *ins[j];
    std::array<std::shared_ptr<const ElementVectorBase>,N_> inDis = ins;
    inDis[j] = stateDis;
    _eval(outRef,inDis);
    VXD difIn(stateDis->getDim());
    VXD difOut(outRef->getDim());
    for(int i=0; i<stateDis->getDim(); i++){
      difIn.setZero();
      difIn(i) = delta;
      ins[j]->boxplus(difIn,stateDis);
      _eval(outDis,inDis);
      outDis->boxminus(outRef,difOut);
      J.col(i) = difOut/delta;
    }
  }

  template<int j>
  bool _testJacInput(const std::array<std::shared_ptr<const ElementVectorBase>,N_>& ins,
      const double& delta = 1e-6, const double& th = 1e-6) const {
    Eigen::MatrixXd J((int)OutPack::d_,(int)std::tuple_element<j,std::tuple<InPacks...>>::type::d_);
    Eigen::MatrixXd J_FD((int)OutPack::d_,(int)std::tuple_element<j,std::tuple<InPacks...>>::type::d_);
    std::shared_ptr<ElementVector> output(new ElementVector(outDefinition_));
    _jac<j>(J,ins);
    _jacFD<j>(J_FD,ins,delta);
    typename Eigen::MatrixXd::Index maxRow, maxCol = 0;
    const double r = (J-J_FD).array().abs().maxCoeff(&maxRow, &maxCol);
    if(r>th){
      std::string outName = outDefinition_->GetName(output->getOuter(maxRow));
      std::string inName = inDefinitions_[j]->GetName(ins[j]->getOuter(maxCol));
      std::cout << "==== Model jacInput (" << j << ") Test failed: " << r
                << " is larger than " << th << " at row "
                << maxRow << "("<< outName << "." << output->getInner(maxRow)
                << ") and col " << maxCol << "("<< inName << "."
                << ins[j]->getInner(maxCol) << ") ====" << std::endl;
      std::cout << "  " << J(maxRow,maxCol) << "  " << J_FD(maxRow,maxCol)
                << std::endl;
      return false;
    } else {
      std::cout << "==== Test successful (" << r << ") ====" << std::endl;
      return true;
    }
  }

  template<int j>
  bool _testJacInput(int& s, const double& delta = 1e-6, const double& th = 1e-6) const{
    std::array<std::shared_ptr<const ElementVectorBase>,N_> ins;
    for(int i=0;i<N_;i++){
      std::shared_ptr<ElementVectorBase> randomState(new ElementVector(inDefinitions_[i]));
      randomState->setRandom(s);
      ins[i] = randomState;
    }
    _testJacInput<j>(ins,delta,th);
  }

  template<int j, int n, int m>
  void _setJacBlock(
      MXD& J, const Eigen::Matrix<double,OutPack::template _GetStateDimension<n>(),
      std::tuple_element<j,std::tuple<InPacks...>>::type::template _GetStateDimension<m>()>& B) const {
    J.block<OutPack::template _GetStateDimension<n>(),
        std::tuple_element<j,std::tuple<InPacks...>>::type::template _GetStateDimension<m>()>(OutPack::template _GetStartIndex<n>(),
            std::tuple_element<j,std::tuple<InPacks...>>::type::template _GetStartIndex<m>()) = B;
  }

 protected:
  std::shared_ptr<ElementVectorDefinition> outDefinition_;
  std::array<std::shared_ptr<ElementVectorDefinition>,N_> inDefinitions_;
};

// ==================== Implementation ==================== //
template<typename ... Ts, typename ... InPacks>
struct TH_pack_size<ElementVectorPack<Ts...>, InPacks...> {
  static constexpr int n_ = TH_pack_size<InPacks...>::n_ + sizeof...(Ts);
};
template<typename ... Ts>
struct TH_pack_size<ElementVectorPack<Ts...>> {
  static constexpr int n_ = sizeof...(Ts);
};

template<int i, typename ... Ts, typename ... Packs>
struct TH_pack_index<i, ElementVectorPack<Ts...>, Packs...> {
  template<int j = i, typename std::enable_if<(j >= sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter() {
    return TH_pack_index<i-sizeof...(Ts),Packs...>::getOuter()+1;
  }
  template<int j=i, typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter() {
    return 0;}
  template<int j=i, typename std::enable_if<(j>=sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner() {
    return TH_pack_index<i-sizeof...(Ts),Packs...>::getInner();}
  template<int j=i, typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner() {
    return i;
  }
};
template<int i, typename ... Ts>
struct TH_pack_index<i, ElementVectorPack<Ts...>> {
  template<typename std::enable_if<(i < sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter() {return 0;};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner() {return i;};
};

}

#endif /* GIF_MODEL_HPP_ */
