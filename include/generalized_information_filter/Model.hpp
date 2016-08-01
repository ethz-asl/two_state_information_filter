/*
 * Model.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_MODEL_HPP_
#define GIF_MODEL_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

class Model{
 public:
  const int n_;
  Model(int n = 1): n_(n){};
  virtual ~Model(){};
  virtual void eval(const State** in, State* out) = 0;
  virtual void jac(const State** in, MXD& out) = 0;
  void jacFD(const State** in, State* out, MXD& J, const double& delta, const int& c = 0){
    assert(c<n_);
    J.setZero();
    const int d = in[c]->d_;
    const State* inDis[n_];
    for(int i=0; i<n_;i++){
       inDis[i] = in[i];
    }
    State* stateDis = in[c]->clone();
    State* outDis = out->clone();
    inDis[c] = stateDis;
    eval(inDis,out);
    VXD difIn(d);
    VXD difOut(out->d_);
    for(int i=0; i<d; i++){
      difIn.setZero();
      difIn(i) = delta;
      in[c]->boxplus(difIn,*stateDis);
      eval(inDis,outDis);
      outDis->boxminus(*out,difOut);
      J.col(i) = difOut/delta;
    }
  }
};

}

#endif /* GIF_MODEL_HPP_ */
