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

  void addRes(std::shared_ptr<BinaryResidualBase> r, const std::string& name = ""){
    binaryResiduals_.push_back(r);
    stateDefinition_->extend(r->preDefinition());
    stateDefinition_->extend(r->posDefinition());
    binaryMeasurements_.push_back(std::list<std::shared_ptr<BinaryMeasurementBase>>());
  }
  void evalRes(const std::shared_ptr<const State>& pre, const std::shared_ptr<const State>& pos){
//    for(int i=0;i<binaryResiduals_.size();i++){
//      auto noi = binaryResiduals_[i]->noiDefinition()->newState();
//      noi->init();
//      auto inn = binaryResiduals_[i]->resDefinition()->newState();
//      binaryResiduals_[i]->evalResidual(inn,preStateWrapper_[i]->wrap(pre),pos,noi);
//      inn->print();
//    }
  }
  std::shared_ptr<StateDefinition> stateDefinition() const{
    return stateDefinition_;
  }

 protected:
  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
  std::vector<std::list<std::shared_ptr<BinaryMeasurementBase>>> binaryMeasurements_;
};

}

#endif /* GIF_FILTER_HPP_ */
