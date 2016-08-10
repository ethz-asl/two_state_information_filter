/*
 * Model.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include <list>
#include <initializer_list>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/StateDefinition.hpp"

namespace GIF{

class TrafoBase{
 public:
  TrafoBase(){};
  virtual ~TrafoBase(){};
  virtual void evalBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut) = 0;
  virtual void jacBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, MXD* J, int i = -1) = 0;
};

template<typename ElementType>
struct DataAndIndex{
  DataAndIndex(const ElementType* x, int i): x_(x), i_(i){};
  const ElementType* x_;
  int i_;
};

template<typename... Ts> class Pack{
 public:
  static constexpr int n_ = sizeof...(Ts);
  typedef std::tuple<Ts...> mtTuple;
  template<int s = 0, int i = 0, size_t n, typename std::enable_if<(i<sizeof...(Ts)-1)>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n>& names, StateDefinition* def){
    typedef typename std::tuple_element<i,std::tuple<Ts...>>::type mtElementType;
    def->addElementDefinition<mtElementType>(names[s+i]);
    addElementToDefinition<s,i+1>(names,def);
  }
  template<int s = 0, int i = 0, size_t n, typename std::enable_if<(i==sizeof...(Ts)-1)>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n>& names, StateDefinition* def){
    typedef typename std::tuple_element<i,std::tuple<Ts...>>::type mtElementType;
    def->addElementDefinition<mtElementType>(names[i]);
  }
};

template<typename... InPacks>
struct TH_pack_size{};

template<typename... Ts, typename... InPacks>
struct TH_pack_size<Pack<Ts...>,InPacks...>{
  static constexpr int n_ = TH_pack_size<InPacks...>::n_+sizeof...(Ts);
};
template<typename... Ts>
struct TH_pack_size<Pack<Ts...>>{
  static constexpr int n_ = sizeof...(Ts);
};

template<int i, typename... InPacks>
struct TH_pack_index{};

template<int i, typename... Ts, typename... InPacks>
struct TH_pack_index<i,Pack<Ts...>,InPacks...>{
  template<typename std::enable_if<(i>=sizeof...(Ts))>::type*>
  static constexpr int getOuter(){return TH_pack_index<i-sizeof...(Ts),InPacks...>::getOuter()+1;};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type*>
  static constexpr int getOuter(){return 0;};
  template<typename std::enable_if<(i>=sizeof...(Ts))>::type*>
  static constexpr int getInner(){return TH_pack_index<i-sizeof...(Ts),InPacks...>::getInner();};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type*>
  static constexpr int getInner(){return i;};
};

template<typename Derived, typename OutPack, typename... InPacks>
class TrafoNew: public TrafoBase{
  static constexpr int n_ = TH_pack_size<InPacks...>::n_;
  static constexpr int m_ = OutPack::n_;

  void evalBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut){
//    _eval(elementBasesIn,elementBasesOut);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    static constexpr int outerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::getOuter();
    static constexpr int innerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::getInner();
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    std::cout << "In: " << std::get<outerIndex>(namesIn_).at(innerIndex) << std::endl;
    _eval(elementBasesIn,elementBasesOut, elements..., dynamic_cast<const Element<mtElementType>*>(elementBasesIn.at(sizeof...(Ps)))->x_);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    static constexpr int innerIndex = sizeof...(Ps)-n_;
    typedef typename OutPack::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    std::cout << "Out: " << namesOut_.at(innerIndex) << std::endl;
    _eval(elementBasesIn,elementBasesOut, elements..., dynamic_cast<Element<mtElementType>*>(elementBasesOut.at(sizeof...(Ps)-n_))->x_);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==n_+m_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    static_cast<Derived&>(*this).eval(elements...);
  }

  const std::tuple<std::array<std::string,InPacks::n_>...> namesIn_;
  const std::array<std::string,m_> namesOut_;
};

template<typename... Packs> class Trafo;

template<typename... Is, typename... Os>
class Trafo<Pack<Is...>, Pack<Os...>>: public TrafoBase{
 public:
  static constexpr int n_ = sizeof...(Is);
  static constexpr int m_ = sizeof...(Os);
  typedef Trafo<Pack<Is...>, Pack<Os...>> mtTrafo;
  Trafo(std::array<std::string,n_> namesIn, std::array<std::string,m_> namesOut): namesIn_(namesIn), namesOut_(namesOut){};
  virtual ~Trafo(){};

  void evalBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut){
    _eval(elementBasesIn,elementBasesOut);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    typedef typename std::tuple_element<sizeof...(Ps),std::tuple<Is...>>::type mtElementType;
    std::cout << "In: " << namesIn_[sizeof...(Ps)] << std::endl;
    _eval(elementBasesIn,elementBasesOut, elements..., dynamic_cast<const Element<mtElementType>*>(elementBasesIn.at(sizeof...(Ps)))->x_);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    typedef typename std::tuple_element<sizeof...(Ps)-n_,std::tuple<Os...>>::type mtElementType;
    std::cout << "Out: " << namesOut_[sizeof...(Ps)-n_] << std::endl;
    _eval(elementBasesIn,elementBasesOut, elements..., dynamic_cast<Element<mtElementType>*>(elementBasesOut.at(sizeof...(Ps)-n_))->x_);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==n_+m_)>::type* = nullptr>
  void _eval(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, Ps&... elements){
    eval(elements...);
  }

  virtual void eval(const Is&... elementsIn, Os&... elementsOut) = 0;

  void jacBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, MXD* J, int i = -1){
    _jacBase(elementBasesIn,elementBasesOut,J,i);
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _jacBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, MXD* J, int i, const DataAndIndex<Ps>&... elements){
    typedef typename std::tuple_element<sizeof...(Ps),std::tuple<Is...>>::type mtElementType;
    _jacBase(elementBasesIn,elementBasesOut,J,i, elements..., DataAndIndex<mtElementType>(&dynamic_cast<const Element<mtElementType>*>(elementBasesIn.at(sizeof...(Ps)))->x_,0));
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _jacBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, MXD* J, int i, const DataAndIndex<Ps>&... elements){
    typedef typename std::tuple_element<sizeof...(Ps)-n_,std::tuple<Os...>>::type mtElementType;
    _jacBase(elementBasesIn,elementBasesOut,J,i, elements..., DataAndIndex<mtElementType>(&dynamic_cast<const Element<mtElementType>*>(elementBasesOut.at(sizeof...(Ps)-n_))->x_,0));
  }

  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==n_+m_)>::type* = nullptr>
  void _jacBase(const std::vector<const ElementBase*>& elementBasesIn, const std::vector<ElementBase*>& elementBasesOut, MXD* J, int i, const DataAndIndex<Ps>&... elements){
    jac(elements...,J,i);
  }

  virtual void jac(const DataAndIndex<Is>&... elementsIn, const DataAndIndex<Os>&... elementsOut, MXD* J, int i = -1) = 0;

  const std::array<std::string,n_> namesIn_;
  const std::array<std::string,m_> namesOut_;
};

template<typename... Is, typename... Js, typename... Packs>
class Trafo<Pack<Is...>, Pack<Js...>, Packs...>: public Trafo<Pack<Is...,Js...>, Packs...>{
 public:
  typedef Trafo<Pack<Is...,Js...>, Packs...> Base;
  typedef Trafo<Pack<Is...>, Pack<Js...>, Packs...> mtTrafo;
  Trafo(std::array<std::string,Base::n_> namesIn, std::array<std::string,Base::m_> namesOut): Base(namesIn,namesOut){};
};

class ResTest: public Trafo<Pack<V3D, QPD>,Pack<V3D>,Pack<V3D>,Pack<V3D>> {
 public:
  ResTest(): mtTrafo({"pos","att","pos","pos"},{"pos"}){};
  void eval(const V3D& posIn, const QPD& attIn, const V3D& posOut, const V3D& posNoi, V3D& posRes){
    std::cout << posIn.transpose() << std::endl;
    std::cout << attIn << std::endl;
  }
  void jac(const DataAndIndex<V3D>& posIn, const DataAndIndex<QPD>& attIn, const DataAndIndex<V3D>& posOut, const DataAndIndex<V3D>& posNoi, const DataAndIndex<V3D>& posRes, MXD* J, int i = -1){
  }
};

class ResTestNew: public TrafoNew<ResTestNew,Pack<V3D>,Pack<V3D, QPD>,Pack<V3D>,Pack<V3D>> {
 public:
  ResTestNew(){};
  void eval(const V3D& posIn, const QPD& attIn, const V3D& posOut, const V3D& posNoi, V3D& posRes){
    std::cout << posIn.transpose() << std::endl;
    std::cout << attIn << std::endl;
  }
//  void jac(const DataAndIndex<V3D>& posIn, const DataAndIndex<QPD>& attIn, const DataAndIndex<V3D>& posOut, const DataAndIndex<V3D>& posNoi, const DataAndIndex<V3D>& posRes, MXD* J, int i = -1){
//  }
};

class Filterr{
 public:
  Filterr(){};
  ~Filterr(){};

  template<typename... Is, typename... Os>
  void addTrafo(Trafo<Pack<Is...>, Pack<Os...>>& t){
    trafos_.push_back(&t);
    Pack<Is...>::addElementToDefinition(t.namesIn_,&inDefinition_);
    Pack<Os...>::addElementToDefinition(t.namesOut_,&outDefinition_);
  }

  template<typename... Pres, typename... Posts, typename... Nois, typename... Ress>
  void addRes(Trafo<Pack<Pres...>, Pack<Posts...>, Pack<Nois...>, Pack<Ress...>>& res){
    res_.push_back(&res);
    Pack<Pres...>::template addElementToDefinition<0>(res.namesIn_,&stateDefinition_);
    Pack<Posts...>::template addElementToDefinition<sizeof...(Pres)>(res.namesIn_,&stateDefinition_);
    Pack<Nois...>::template addElementToDefinition<sizeof...(Pres)+sizeof...(Posts)>(res.namesIn_,&noiseDefinition_);
    Pack<Ress...>::addElementToDefinition(res.namesOut_,&resDefinition_);
  }

  void eval(const State* in,State* out){
    for(auto t : trafos_){
      t->evalBase(in->getElements(), out->getElements());
    }
  }

  void evalResidual(const State* pre,const State* post,const State* noi,State* res){
    std::vector<const ElementBase*> in;
    std::vector<const ElementBase*> preVec = pre->getElements();
    std::vector<const ElementBase*> postVec = post->getElements();
    std::vector<const ElementBase*> noiVec = noi->getElements();
    in.insert(in.end(),preVec.begin(),preVec.end());
    in.insert(in.end(),postVec.begin(),postVec.end());
    in.insert(in.end(),noiVec.begin(),noiVec.end());
    for(auto t : res_){
      t->evalBase(in, res->getElements());
    }
  }
  StateDefinition inDefinition_;
  StateDefinition outDefinition_;
  StateDefinition stateDefinition_;
  StateDefinition noiseDefinition_;
  StateDefinition resDefinition_;
  std::vector<TrafoBase*> trafos_;
  std::vector<TrafoBase*> res_;
};

class TrafoTest: public Trafo<Pack<V3D, QPD>,Pack<V3D>> {
 public:
  typedef Trafo<Pack<V3D, QPD>,Pack<V3D>> Base;
  TrafoTest(std::array<std::string,n_> namesIn, std::array<std::string,m_> namesOut): Base(namesIn,namesOut){};
  void eval(const V3D& pos, const QPD& att, V3D& posOut){
    std::cout << pos.transpose() << std::endl;
    std::cout << att << std::endl;
  }
  void jac(const DataAndIndex<V3D>& pos, const DataAndIndex<QPD>& att, const DataAndIndex<V3D>& posOut, MXD* J, int i = -1){
  }
};

class Model{
 public:
  Model(){};
  virtual ~Model(){};
  virtual void _eval(const std::vector<const State*>& in, State* out) = 0;
  virtual void _jac(const std::vector<const State*>& in, MXD& out, int c) = 0;
  void _jacFD(const std::vector<const State*>& in, MXD& J, int c, const double& delta){
    J.setZero();
    State* stateDis = inputDefinitions_[c]->newState();
    *stateDis = *in[c];
    std::vector<const State*> inDis(in);
    inDis[c] = stateDis;
    State* outRef = outputDefinition_->newState();
    State* outDis = outputDefinition_->newState();
    _eval(inDis,outRef);
    VXD difIn(inputDefinitions_[c]->getDim());
    VXD difOut(outputDefinition_->getDim());
    for(int i=0; i<inputDefinitions_[c]->getDim(); i++){
      difIn.setZero();
      difIn(i) = delta;
      inputDefinitions_[c]->boxplus(in[c],difIn,stateDis);
      _eval(inDis,outDis);
      outputDefinition_->boxminus(outDis,outRef,difOut);
      J.col(i) = difOut/delta;
    }
    delete stateDis;
    delete outRef;
    delete outDis;
  }
  void initStateDefinitions(std::shared_ptr<StateDefinition> outputDefinition, std::initializer_list<std::shared_ptr<StateDefinition>> inputList){
    outputDefinition_ = outputDefinition;
    inputDefinitions_ = inputList;
    buildStateDefinitions();
  }

 protected:
  virtual void buildStateDefinitions() = 0;

 private:
  std::vector<std::shared_ptr<StateDefinition>> inputDefinitions_;
  std::shared_ptr<StateDefinition> outputDefinition_;
};

}

#endif /* GIF_MODEL_HPP_ */
