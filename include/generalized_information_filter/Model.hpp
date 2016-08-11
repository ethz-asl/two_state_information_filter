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

template<typename... Ts>
class Pack{
 public:
  static constexpr int n_ = sizeof...(Ts);
  typedef std::tuple<Ts...> mtTuple;
  template<int i = 0, typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n_>& names, std::shared_ptr<StateDefinition> def){
    typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
    def->addElementDefinition<mtElementType>(names.at(i));
    addElementToDefinition<i+1>(names,def);
  }
  template<int i = 0, typename std::enable_if<(i==sizeof...(Ts))>::type* = nullptr>
  static void addElementToDefinition(const std::array<std::string,n_>& names, std::shared_ptr<StateDefinition> def){}

  template<int i = 0, typename std::enable_if<(i<sizeof...(Ts))>::type* = nullptr>
  static void addElementToDefinition(std::shared_ptr<StateDefinition> def){
    typedef typename std::tuple_element<i,mtTuple>::type mtElementType;
    def->addElementDefinition<mtElementType>(std::to_string(i));
    addElementToDefinition<i+1>(def);
  }
  template<int i = 0, typename std::enable_if<(i==sizeof...(Ts))>::type* = nullptr>
  static void addElementToDefinition(std::shared_ptr<StateDefinition> def){}
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

class ModelBase{
 public:
  ModelBase(){};
  virtual ~ModelBase(){};
  virtual void evalBase(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn) = 0;
};

template<typename Derived, typename OutPack, typename... InPacks>
class Model: public ModelBase{
 public:
  static constexpr int n_ = OutPack::n_;
  static constexpr int m_ = TH_pack_size<InPacks...>::n_;
  static constexpr int N_ = sizeof...(InPacks);
  typedef Model<Derived,OutPack,InPacks...> mtBase;
  Model(){}

  void evalBase(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn){
    _eval(elementBasesOut,elementBasesIn);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)<n_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static constexpr int innerIndex = sizeof...(Ps);
    typedef typename OutPack::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _eval(elementBasesOut,elementBasesIn, elements..., dynamic_cast<Element<mtElementType>*>(elementBasesOut.at(sizeof...(Ps)))->x_);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)>=n_ & sizeof...(Ps)<n_+m_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static constexpr int outerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getOuter();
    static constexpr int innerIndex = TH_pack_index<sizeof...(Ps)-n_,InPacks...>::getInner();
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _eval(elementBasesOut,elementBasesIn, elements..., dynamic_cast<const Element<mtElementType>*>(elementBasesIn.at(sizeof...(Ps)-n_))->x_);
  }
  template<typename... Ps, typename std::enable_if<(sizeof...(Ps)==m_+n_)>::type* = nullptr>
  void _eval(const std::vector<ElementBase*>& elementBasesOut, const std::vector<const ElementBase*>& elementBasesIn, Ps&... elements){
    static_cast<Derived&>(*this).eval(elements...);
  }

  template<int n, typename... Ps, typename std::enable_if<(sizeof...(Ps)<m_)>::type* = nullptr>
  void _jac(MXD* J, const std::vector<std::pair<const ElementBase*,int>>& elementsIn, Ps&... elements){
    static constexpr int outerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::getOuter();
    static constexpr int innerIndex = TH_pack_index<sizeof...(Ps),InPacks...>::getInner();
    typedef typename std::tuple_element<outerIndex,std::tuple<InPacks...>>::type::mtTuple mtTuple;
    typedef typename std::tuple_element<innerIndex,mtTuple>::type mtElementType;
    _jac(J,elementsIn, elements..., dynamic_cast<const Element<mtElementType>*>(elementsIn.at(sizeof...(Ps)).first)->x_);
  }
  template<int n, typename... Ps, typename std::enable_if<(sizeof...(Ps)==m_)>::type* = nullptr>
  void _jac(MXD* J, const std::vector<std::pair<const ElementBase*,int>>& elementsIn, Ps&... elements){
    static_assert(n<N_,"No such Jacobian!");
    static_cast<Derived&>(*this).template eval<n>(J,elements...);
  }
};

}

#endif /* GIF_MODEL_HPP_ */
