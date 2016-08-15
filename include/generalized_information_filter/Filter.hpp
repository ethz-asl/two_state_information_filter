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

class Filter{
 public:
  Filter(): stateDefinition_(new StateDefinition()){};
  virtual ~Filter(){};

  void addRes(std::shared_ptr<BinaryResidualBase> r){
    binaryResiduals_.push_back(r);
    stateDefinition_->extend(r->preDefinition());
    stateDefinition_->extend(r->posDefinition());
  }

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

  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
};

}

#endif /* GIF_FILTER_HPP_ */
