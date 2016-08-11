#include "gtest/gtest.h"

#include "generalized_information_filter/Filter.hpp"

using namespace GIF;

//class OdometryResidual: public PairwiseResidual{
// public:
//  OdometryResidual(){};
//  virtual ~OdometryResidual(){};
//  void eval(const State* inPre, const State* inPost, const State* noise, State* res){
////    Eigen::Matrix<double,5,1> v;
////    v << 1,2,3,4,5;
////    posOut_->get(out) = (posIn_->get(in)+v)*(timeOffset_->get(in)+2.0);
//  }
//  void jacPre(const State* inPre, const State* inPost, MXD& out){
////    Eigen::Matrix<double,5,1> v;
////    v << 1,2,3,4,5;
////    out.setZero();
////    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()).setIdentity();
////    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()) *= (timeOffset_->get(in)+2.0);
////    out.block(posOut_->getIndex(),timeOffset_->getIndex(),posOut_->getDim(),timeOffset_->getDim()) = posIn_->get(in)+v;
//  }
//  void jacPost(const State* inPre, const State* inPost, MXD& out){
//
//  }
//  void jacNoise(const State* inPre, const State* inPost, MXD& out){
//
//  }
// protected:
//  ElementDefinition<double>* yawState_ = nullptr;
//  ElementDefinition<Eigen::Matrix<double,2,1>>* posState_ = nullptr;
//  ElementDefinition<double>* yawRes_ = nullptr;
//  ElementDefinition<Eigen::Matrix<double,2,1>>* posRes_ = nullptr;
//  ElementDefinition<double>* yawNoise_ = nullptr;
//  ElementDefinition<Eigen::Matrix<double,2,1>>* posNoise_ = nullptr;
//  ElementDefinition<double>* rorMeas_ = nullptr;
//  ElementDefinition<Eigen::Matrix<double,2,1>>* velMeas_ = nullptr;
//
//  virtual void buildStateDefinitions(){
//    yawState_ = stateDefinition_->addElementDefinition<double>("yaw");
//    posState_ = stateDefinition_->addElementDefinition<Eigen::Matrix<double,2,1>>("pos");
//    yawRes_ = residualDefinition_->addElementDefinition<double>("yaw");
//    posRes_ = residualDefinition_->addElementDefinition<Eigen::Matrix<double,2,1>>("pos");
//    yawNoise_ = noiseDefinition_->addElementDefinition<double>("yaw");
//    posNoise_ = noiseDefinition_->addElementDefinition<Eigen::Matrix<double,2,1>>("pos");
//    rorMeas_ = measDefinition_->addElementDefinition<double>("ror");
//    velMeas_ = measDefinition_->addElementDefinition<Eigen::Matrix<double,2,1>>("vel");
//  }
//};

// The fixture for testing class ScalarState
class FilterTest : public virtual ::testing::Test {
 protected:
  FilterTest(){}
  virtual ~FilterTest(){}
};

// Test constructors
TEST_F(FilterTest, constructor) {
//  std::cout << "================================" << std::endl;
//  Filter f;
//  ResTest res({"pos"},{"pos","att","pos","pos"});
//  f.addRes(res);
//  State* pre = f.stateDefinition_.newState();
//  State* post = f.stateDefinition_.newState();
//  State* noi = f.noiseDefinition_.newState();
//  State* resi = f.resDefinition_.newState();
//  f.evalResidual(pre,post,noi,resi);
//  std::cout << "================================" << std::endl;




//  // Odometry residuals
//  OdometryResidual odometryResidual1;
//  OdometryResidual odometryResidual2;
//
//  Filter filter;
//  int odometryResidual1Ind = filter.addResidual(&odometryResidual1);
//  int odometryResidual2Ind = filter.addResidual(&odometryResidual2);
//  filter.init();
//
//  filter.stateDefinition()->print(filter.getState());

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
