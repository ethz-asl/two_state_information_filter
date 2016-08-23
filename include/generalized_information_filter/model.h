#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include <list>
#include <initializer_list>

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/state-definition.h"

namespace GIF {

template<typename ... Ts>
struct TH_pack_dim;
template<typename T, typename ... Ts>
struct TH_pack_dim<T, Ts...> {
  static constexpr int d_ = TH_pack_dim<Ts...>::d_ + ElementTraits<T>::d_;
};
template<typename T>
struct TH_pack_dim<T> {
  static constexpr int d_ = ElementTraits<T>::d_;
};
template<>
struct TH_pack_dim<> {
  static constexpr int d_ = 0;
};

template<typename ... Ts>
class ElementPack {
 public:
  static constexpr int n_ = sizeof...(Ts);
  static constexpr int d_ = TH_pack_dim<Ts...>::d_;
  typedef std::tuple<Ts...> mtTuple;

  static std::shared_ptr<StateDefinition> makeStateDefinition(const std::array<std::string,n_>& names) {
    std::shared_ptr<StateDefinition> def(new StateDefinition());
    addElementsToDefinition(names,def);
    return def;
  }

  template<int i = 0, typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static void addElementsToDefinition(const std::array<std::string,n_>& names,
      const std::shared_ptr<StateDefinition>& def) {
    typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
    def->addElementDefinition<mtElementType>(names.at(i));
    addElementsToDefinition<i+1>(names,def);
  }
  template<int i = 0, typename std::enable_if<(i==sizeof...(Ts))>::type* = nullptr>
  static void addElementsToDefinition(const std::array<std::string,n_>& names,
      const std::shared_ptr<StateDefinition>& def) {}

  template<int i>
  static constexpr int getDim() {
    return ElementTraits<typename std::tuple_element<i,mtTuple>::type>::d_;
  }

  template<int i, typename std::enable_if<(i>0)>::type* = nullptr>
  static constexpr int getStart() {
    return getStart<i-1>() + getDim<i-1>();
  }
  template<int i = 0, typename std::enable_if<(i==0)>::type* = nullptr>
  static constexpr int getStart() {
    return 0;
  }
};

template<typename ... InPacks>
struct TH_pack_size;
template<typename ... Ts, typename ... InPacks>
struct TH_pack_size<ElementPack<Ts...>, InPacks...> {
  static constexpr int n_ = TH_pack_size<InPacks...>::n_ + sizeof...(Ts);
};
template<typename ... Ts>
struct TH_pack_size<ElementPack<Ts...>> {
  static constexpr int n_ = sizeof...(Ts);
};

template<int i, typename ... Packs>
struct TH_pack_index;
template<int i, typename ... Ts, typename ... Packs>
struct TH_pack_index<i, ElementPack<Ts...>, Packs...> {
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
struct TH_pack_index<i, ElementPack<Ts...>> {
  template<typename std::enable_if<(i < sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter() {return 0;};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner() {return i;};
};

template<typename Derived, typename OutPack, typename ... InPacks>
class Model {
 public:
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  static constexpr int N_ = sizeof...(InPacks);
  typedef Model<Derived,OutPack,InPacks...> mtBase;

  Model(const std::array<std::string,n_>& namesOut,
        const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    outDefinition_ = OutPack::makeStateDefinition(namesOut);
    makeInDefinitons(namesIn);
  }

  template<int i = 0, typename std::enable_if<(i<N_)>::type* = nullptr>
  void makeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {
    inDefinitions_[i] =
        std::tuple_element<i,std::tuple<InPacks...>>::type::makeStateDefinition(
            std::get<i>(namesIn));
    makeInDefinitons<i+1>(namesIn);
  }

  template<int i = 0, typename std::enable_if<(i==N_)>::type* = nullptr>
  void makeInDefinitons(
      const std::tuple<std::array<std::string,InPacks::n_>...>& namesIn) {}

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::shared_ptr<StateBase>& out,
             const std::array<std::shared_ptr<const StateBase>,N_>& ins,
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
  void _eval(const std::shared_ptr<StateBase>& out,
             const std::array<std::shared_ptr<const StateBase>,N_>& ins,
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
  void _eval(const std::shared_ptr<StateBase>& out,
             const std::array<std::shared_ptr<const StateBase>,N_>& ins,
             Ps&... elements) const{
    static_cast<const Derived&>(*this).eval(elements...);
  }

  template<int j, typename... Ps, typename std::enable_if<(sizeof...(Ps)<m_)>::type* = nullptr>
  void _jac(MXD& J, const std::array<std::shared_ptr<const StateBase>,N_>& ins,
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
  void _jac(MXD& J, const std::array<std::shared_ptr<const StateBase>,N_>& ins,
            Ps&... elements) const{
    static_assert(j<N_,"No such Jacobian!");
    static_cast<const Derived&>(*this).template jac<j>(J,elements...);
  }

  template<int j>
  void _jacFD(MXD& J,
              const std::array<std::shared_ptr<const StateBase>,N_>& ins,
              const double& delta = 1e-8) const{
    std::shared_ptr<State> stateDis(new State(inDefinitions_[j]));
    std::shared_ptr<State> outRef(new State(outDefinition_));
    std::shared_ptr<State> outDis(new State(outDefinition_));
    J.resize(outRef->getDim(),stateDis->getDim());
    J.setZero();
    *stateDis = *ins[j];
    std::array<std::shared_ptr<const StateBase>,N_> inDis = ins;
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
  bool _testJacInput(const std::array<std::shared_ptr<const StateBase>,N_>& ins,
      const double& delta = 1e-6, const double& th = 1e-6) const {
    Eigen::MatrixXd J((int)OutPack::d_,(int)std::tuple_element<j,std::tuple<InPacks...>>::type::d_);
    Eigen::MatrixXd J_FD((int)OutPack::d_,(int)std::tuple_element<j,std::tuple<InPacks...>>::type::d_);
    std::shared_ptr<State> output(new State(outDefinition_));
    _jac<j>(J,ins);
    _jacFD<j>(J_FD,ins,delta);
    typename Eigen::MatrixXd::Index maxRow, maxCol = 0;
    const double r = (J-J_FD).array().abs().maxCoeff(&maxRow, &maxCol);
    if(r>th){
      std::string outName = outDefinition_->getName(output->getOuter(maxRow));
      std::string inName = inDefinitions_[j]->getName(ins[j]->getOuter(maxCol));
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
    std::array<std::shared_ptr<const StateBase>,N_> ins;
    for(int i=0;i<N_;i++){
      std::shared_ptr<StateBase> randomState(new State(inDefinitions_[i]));
      randomState->setRandom(s);
      ins[i] = randomState;
    }
    _testJacInput<j>(ins,delta,th);
  }

  template<int j, int n, int m>
  void _setJacBlock(
      MXD& J, const Eigen::Matrix<double,OutPack::template getDim<n>(),
      std::tuple_element<j,std::tuple<InPacks...>>::type::template getDim<m>()>& B) const {
    J.block<OutPack::template getDim<n>(),
        std::tuple_element<j,std::tuple<InPacks...>>::type::template getDim<m>()>(OutPack::template getStart<n>(),
            std::tuple_element<j,std::tuple<InPacks...>>::type::template getStart<m>()) = B;
  }

 protected:
  std::shared_ptr<StateDefinition> outDefinition_;
  std::array<std::shared_ptr<StateDefinition>,N_> inDefinitions_;
};

}

#endif /* GIF_MODEL_HPP_ */
