#include "gtest/gtest.h"
#include <assert.h>

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/Transformation.hpp"

using namespace GIF;

class State1: public State{
 public:
  State1(): State(){
    timeOffset_ = addElement<double>();
    pos_ = addElement<Eigen::Matrix<double,5,1>>();
    fea_ = addElement<std::array<Eigen::Matrix<double,3,1>,3>>();
  }
  State* clone() const{
    State1* s = new State1();
    *dynamic_cast<State*>(s) = *dynamic_cast<const State*>(this);
    return s;
  }
  Element<double>* timeOffset_;
  Element<Eigen::Matrix<double,5,1>>* pos_;
  Element<std::array<Eigen::Matrix<double,3,1>,3>>* fea_;
};

class TransformationExample: public Transformation{
 public:
  TransformationExample(){};
  virtual ~TransformationExample(){};
  void eval(const State** in, State* out){
    const State1* inCast = dynamic_cast<const State1*>(*in);
    State1* outCast = dynamic_cast<State1*>(out);
    outCast->pos_->x_ = inCast->pos_->x_;
  }
  void jac(const State** in, MXD& out){
    const State1* inCast = dynamic_cast<const State1*>(*in);
    out.setZero();
    out.block(inCast->pos_->i_,inCast->pos_->i_,inCast->pos_->d_,inCast->pos_->d_).setIdentity();
  }
};

// The fixture for testing class ScalarState
class NewStateTest : public virtual ::testing::Test {
 protected:
  NewStateTest():covMat_(1,1) {
//    unsigned int s = 1;
//    testElement1_.setRandom(s);
//    testElement2_.setRandom(s);
  }
  virtual ~NewStateTest() {
  }
//  GIF::ScalarElement testElement1_;
//  GIF::ScalarElement testElement2_;
//  GIF::ScalarElement::mtDifVec difVec_;
  MXD covMat_;
};

// Test constructors
TEST_F(NewStateTest, constructor) {
  State1 s, s2;
  s.print();
  Eigen::VectorXd v(s.d_);
  v.setZero();
  for(int i = 0; i < s.d_; i++){
    v(i) = i;
  }
  s.boxplus(v,s2);
  s2.print();
  s.boxminus(s2,v);
  std::cout << v.transpose() << std::endl;


  TransformationExample t;
  const State* states[1];
  states[0] = &s;
  MXD J(s.d_,s.d_);
  t.jac(states,J);
  std::cout << J << std::endl;
  MXD JFD(s.d_,s.d_);
  t.jacFD(states,&s2,JFD,1e-8);
  std::cout << JFD << std::endl;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
