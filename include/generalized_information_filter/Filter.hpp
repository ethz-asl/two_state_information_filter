/*
 * Filter.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include <map>

#include "generalized_information_filter/BinaryResidual.hpp"
#include "generalized_information_filter/common.hpp"

namespace GIF{

//class Filter{
// public:
//  Filter(){};
//  ~Filter(){};
//
//  template<typename T, typename... Ress, typename... Pres, typename... Posts, typename... Nois>
//  void addRes(Model<T, Pack<Ress...>, Pack<Pres...>, Pack<Posts...>, Pack<Nois...>>& res){
//    res_.push_back(&res);
//    Pack<Ress...>::addElementToDefinition(res.namesOut_,&resDefinition_);
//    Pack<Pres...>::addElementToDefinition(std::get<0>(res.namesIn_),&stateDefinition_);
//    Pack<Posts...>::addElementToDefinition(std::get<1>(res.namesIn_),&stateDefinition_);
//    Pack<Nois...>::addElementToDefinition(std::get<2>(res.namesIn_),&noiseDefinition_);
//  }
//
//  void evalResidual(const State* pre,const State* post,const State* noi,State* res){
//    std::vector<const ElementBase*> in;
//    std::vector<const ElementBase*> preVec = pre->getElements();
//    std::vector<const ElementBase*> postVec = post->getElements();
//    std::vector<const ElementBase*> noiVec = noi->getElements();
//    in.insert(in.end(),preVec.begin(),preVec.end());
//    in.insert(in.end(),postVec.begin(),postVec.end());
//    in.insert(in.end(),noiVec.begin(),noiVec.end());
//    for(auto t : res_){
//      t->evalBase(res->getElements(), in);
//    }
//  }
//  StateDefinition stateDefinition_;
//  StateDefinition noiseDefinition_;
//  StateDefinition resDefinition_;
//  std::vector<ModelBase*> res_;
//};

}

#endif /* GIF_FILTER_HPP_ */
