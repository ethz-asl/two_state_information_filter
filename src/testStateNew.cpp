#include "gtest/gtest.h"
#include <assert.h>
#include <array>

#include "../include/generalized_information_filter/PairwiseResidual.hpp"
#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/Transformation.hpp"

using namespace GIF;

class TransformationExample: public Transformation{
 public:
  TransformationExample(StateDefinition* inputDefinition, StateDefinition* outputDefinition):
    Transformation(inputDefinition,outputDefinition){
    timeOffset_ = inputDefinition->addElementDefinition<double>("timeOffset");
    posIn_ = inputDefinition->addElementDefinition<Eigen::Matrix<double,5,1>>("pos");
    fea_ = inputDefinition->addElementDefinition<std::array<Eigen::Matrix<double,3,1>,3>>("fea");
    posOut_ = outputDefinition->addElementDefinition<Eigen::Matrix<double,5,1>>("pos");
  };
  virtual ~TransformationExample(){};
  void eval(const State* in, State* out){
    Eigen::Matrix<double,5,1> v;
    v << 1,2,3,4,5;
    posOut_->get(out) = (posIn_->get(in)+v)*(timeOffset_->get(in)+2.0);
  }
  void jac(const State* in, MXD& out){
    Eigen::Matrix<double,5,1> v;
    v << 1,2,3,4,5;
    out.setZero();
    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()).setIdentity();
    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()) *= (timeOffset_->get(in)+2.0);
    out.block(posOut_->getIndex(),timeOffset_->getIndex(),posOut_->getDim(),timeOffset_->getDim()) = posIn_->get(in)+v;
  }
 private:
  ElementDefinition<double>* timeOffset_;
  ElementDefinition<Eigen::Matrix<double,5,1>>* posIn_;
  ElementDefinition<std::array<Eigen::Matrix<double,3,1>,3>>* fea_;
  ElementDefinition<Eigen::Matrix<double,5,1>>* posOut_;
};

// The fixture for testing class ScalarState
class NewStateTest : public virtual ::testing::Test {
 protected:
  NewStateTest():covMat_(1,1) {
  }
  virtual ~NewStateTest() {
  }
  MXD covMat_;
};

// Test constructors
TEST_F(NewStateTest, constructor) {
  StateDefinition def1;
  StateDefinition def2;
  TransformationExample t(&def1,&def2);
  State* s1a = def1.newState();
  State* s1b = def1.newState();
  def1.init(s1a);
  def1.print(s1a);

  // Boxplus and boxminus
  Eigen::VectorXd v(def1.getDim());
  v.setZero();
  for(int i = 0; i < def1.getDim(); i++){
    v(i) = i;
  }
  def1.boxplus(s1a,v,s1b);
  def1.print(s1b);
  def1.boxminus(s1a,s1b,v);
  std::cout << v.transpose() << std::endl;

  // Transformation
  State* s2 = def2.newState();
  MXD J(def2.getDim(),def1.getDim());
  t.jac(s1a,J);
  std::cout << J << std::endl;
  MXD JFD(def2.getDim(),def1.getDim());
  t.jacFD(s1a,JFD,1e-8);
  std::cout << JFD << std::endl;

  delete s1a;
  delete s1b;
//  delete s2;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
