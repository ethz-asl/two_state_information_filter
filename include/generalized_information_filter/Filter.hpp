/*
 * Filter.hpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_FILTER_HPP_
#define GIF_FILTER_HPP_

#include <map>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/PairwiseResidual.hpp"

namespace GIF{

class MeasurementList{
 public:
  MeasurementList(){};
  void addMeasurement(const State* meas, const double& time){
    std::pair<std::map<double,const State*>::iterator,bool> ret;
    ret = measurements_.insert(std::pair<double,const State*>(time,meas));
    if (ret.second==false) {
      std::cout << "ERROR: measurement already existed";
    }
  }

 protected:
  std::map<double, const State*> measurements_;
};

class Filter{
 public:
  Filter(): state_(nullptr), time_(0), cov_(0,0), initializedTime_(false), stateDefinition_(new StateDefinition()){};
  ~Filter(){};
  int addResidual(PairwiseResidual* res){
    residuals_.push_back(std::pair<PairwiseResidual*, MeasurementList>(res,MeasurementList()));
    return residuals_.size()-1;
  }
  void addMeasurement(int i, const State* meas, const double& time){
    residuals_[i].second.addMeasurement(meas,time);
  }
  void init(){
    for(auto r : residuals_){
      r.first->initStateDefinitions(stateDefinition_);
    }
    time_ = 0;
    initializedTime_ = false;
    state_ = stateDefinition_->newState();
    cov_.resize(stateDefinition_->getDim(),stateDefinition_->getDim());
  }
  std::shared_ptr<const StateDefinition> stateDefinition(){
    return stateDefinition_;
  }
  void updateStep(const double& time){
    // check involved residuals

    // construct residual and Jacobian

    // do GIF update

  }
  void update(){
    // compute maximal update time

    // split and merge (different modes)

    // execute single update steps

  }
  const State* getState(){
    return state_;
  }

 protected:
  std::vector<std::pair<PairwiseResidual*, MeasurementList>> residuals_;
  State* state_;
  double time_;
  MXD cov_;
  bool initializedTime_;
  std::shared_ptr<StateDefinition> stateDefinition_;

};

}

#endif /* GIF_FILTER_HPP_ */
