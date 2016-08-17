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

  void addRes(const std::shared_ptr<BinaryResidualBase>& r, const std::string& name = ""){
    binaryResiduals_.push_back(r);
    stateDefinition_->extend(r->preDefinition());
    stateDefinition_->extend(r->posDefinition());
    binaryMeasurements_.push_back(std::list<std::shared_ptr<MeasurementBase>>());
    binaryWrappersPre_.emplace_back(new StateWrapper(r->preDefinition(),stateDefinition_));
    binaryWrappersPos_.emplace_back(new StateWrapper(r->posDefinition(),stateDefinition_));
  }
  void evalRes(const std::shared_ptr<const StateBase>& pre, const std::shared_ptr<const StateBase>& pos){
    for(int i=0;i<binaryResiduals_.size();i++){
      std::shared_ptr<State> inn(new State(binaryResiduals_[i]->resDefinition()));
      std::shared_ptr<State> noi(new State(binaryResiduals_[i]->noiDefinition()));
      noi->setIdentity();
      binaryWrappersPre_[i]->setState(pre);
      binaryWrappersPos_[i]->setState(pos);
      binaryResiduals_[i]->evalResidual(inn,binaryWrappersPre_[i],binaryWrappersPos_[i],noi);
      inn->print();
    }
  }
  std::shared_ptr<StateDefinition> stateDefinition() const{
    return stateDefinition_;
  }

 protected:
  std::shared_ptr<StateDefinition> stateDefinition_;
  std::vector<std::shared_ptr<BinaryResidualBase>> binaryResiduals_;
  std::vector<std::list<std::shared_ptr<MeasurementBase>>> binaryMeasurements_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPre_;
  std::vector<std::shared_ptr<const StateWrapper>> binaryWrappersPos_;
};

}

#endif /* GIF_FILTER_HPP_ */
