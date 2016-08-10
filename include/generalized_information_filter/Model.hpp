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
  virtual void evalBase(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn) = 0;
//  virtual void jacBase(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, MXD* J, int i = -1) = 0;
};

//template<typename ElementType>
//struct DataAndIndex{
//  DataAndIndex(const ElementType* x, int i): x_(x), i_(i){};
//  const ElementType* x_;
//  int i_;
//};

template<typename... Ts> class Pack{
 public:
  static constexpr int n_ = sizeof...(Ts);
  typedef std::tuple<Ts...> mtTuple;
  template<int i = 0, typename std::enable_if<(i<sizeof...(Ts)-1)>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n_>& names, StateDefinition* def){
    typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
    def->addElementDefinition<mtElementType>(names.at(i));
    addElementToDefinition<i+1>(names,def);
  }
  template<int i = 0, typename std::enable_if<(i==sizeof...(Ts)-1)>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n_>& names, StateDefinition* def){
    typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
    def->addElementDefinition<mtElementType>(names.at(i));
  }
};

template<typename... InPacks>
struct TH_pack_size;
template<typename... Ts, typename... InPacks>
struct TH_pack_size<Pack<Ts...>,InPacks...>{
  static constexpr int n_ = TH_pack_size<InPacks...>::n_+sizeof...(Ts);
};
template<typename... Ts>
struct TH_pack_size<Pack<Ts...>>{
  static constexpr int n_ = sizeof...(Ts);
};

template<int i, typename... Packs>
struct TH_pack_index;
template<int i, typename... Ts, typename... Packs>
struct TH_pack_index<i,Pack<Ts...>,Packs...>{
  template<int j=i, typename std::enable_if<(j>=sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter(){return TH_pack_index<i-sizeof...(Ts),Packs...>::getOuter()+1;}
  template<int j=i, typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter(){return 0;}
  template<int j=i, typename std::enable_if<(j>=sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner(){return TH_pack_index<i-sizeof...(Ts),Packs...>::getInner();}
  template<int j=i, typename std::enable_if<(j<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner(){return i;}
};
template<int i, typename... Ts>
struct TH_pack_index<i,Pack<Ts...>>{
  template<typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getOuter(){return 0;};
  template<typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static constexpr int getInner(){return i;};
};

template<typename Derived, typename OutPack, typename... InPacks>
class Trafo: public TrafoBase{
 public:
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  typedef Trafo<Derived,OutPack,InPacks...> mtTrafo;
  Trafo(std::array<std::string,n_> namesOut, std::array<std::string,m_> namesIn): namesOut_(namesOut){
    _initNamesIn(namesIn);
  };

  template<int i=0, typename std::enable_if<(i<m_)>::type* = nullptr>
  void _initNamesIn(std::array<std::string,m_> namesIn){
    static constexpr int outerIndex = TH_pack_index<i,InPacks...>::getOuter();
    static constexpr int innerIndex = TH_pack_index<i,InPacks...>::getInner();
    std::get<outerIndex>(namesIn_).at(innerIndex) = namesIn.at(i);
    _initNamesIn<i+1>(namesIn);
  }
  template<int i=0, typename std::enable_if<(i==m_)>::type* = nullptr>
  void _initNamesIn(std::array<std::string,m_> namesIn){}

  void evalBase(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn){
    _eval(elementBasesOut,elementBasesIn);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static constexpr int innerIndex = sizeof...(Ps);
    typedef typename OutPack::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    std::cout << "Out: " << namesOut_.at(innerIndex) << std::endl;
    _eval(elementBasesOut,elementBasesIn, elements..., dynamic_cast<Element<mtElementType>*>(elementBasesOut.at(sizeof...(Ps)))->x_);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static constexpr int outerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getOuter();
    static constexpr int innerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getInner();
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    std::cout << "In: " << std::get<outerIndex>(namesIn_).at(innerIndex) << std::endl;
    _eval(elementBasesOut,elementBasesIn, elements..., dynamic_cast<const Element<mtElementType>*>(elementBasesIn.at(sizeof...(Ps)-n_))->x_);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==m_+n_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static_cast<Derived&>(*this).eval(elements...);
  }

  std::array<std::string,n_> namesOut_;
  std::tuple<std::array<std::string,InPacks::n_>...> namesIn_;
};

class TrafoTest: public Trafo<TrafoTest,Pack<V3D>,Pack<V3D, QPD>> {
 public:
  TrafoTest(std::array<std::string,n_> namesOut, std::array<std::string,m_> namesIn): mtTrafo(namesOut,namesIn){};
  void eval(V3D& posOut, const V3D& pos, const QPD& att){
    std::cout << pos.transpose() << std::endl;
    std::cout << att << std::endl;
  }
//  void jac(const DataAndIndex<V3D>& pos, const DataAndIndex<QPD>& att, const DataAndIndex<V3D>& posOut, MXD* J, int i = -1){
//  }
};

class ResTest: public Trafo<ResTest,Pack<V3D>,Pack<V3D, QPD>,Pack<V3D>,Pack<V3D>> {
 public:
  ResTest(std::array<std::string,n_> namesOut, std::array<std::string,m_> namesIn): mtTrafo(namesOut,namesIn){};
  void eval(V3D& posRes, const V3D& posIn, const QPD& attIn, const V3D& posOut, const V3D& posNoi){
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

  template<typename T, typename... Ress, typename... Pres, typename... Posts, typename... Nois>
  void addRes(Trafo<T, Pack<Ress...>, Pack<Pres...>, Pack<Posts...>, Pack<Nois...>>& res){
    res_.push_back(&res);
    Pack<Ress...>::addElementToDefinition(res.namesOut_,&resDefinition_);
    Pack<Pres...>::addElementToDefinition(std::get<0>(res.namesIn_),&stateDefinition_);
    Pack<Posts...>::addElementToDefinition(std::get<1>(res.namesIn_),&stateDefinition_);
    Pack<Nois...>::addElementToDefinition(std::get<2>(res.namesIn_),&noiseDefinition_);
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
      t->evalBase(res->getElements(), in);
    }
  }
  StateDefinition stateDefinition_;
  StateDefinition noiseDefinition_;
  StateDefinition resDefinition_;
  std::vector<TrafoBase*> res_;
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
