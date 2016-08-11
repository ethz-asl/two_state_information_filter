#include "gtest/gtest.h"

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/Transformation.hpp"

using namespace GIF;

class TransformationExample: public Transformation<ElementPack<V3D>,ElementPack<double, std::array<V3D,4>>>{
 public:
  TransformationExample(){};
  virtual ~TransformationExample(){};
  void evalTransform(V3D& posOut,const double& timeIn,const std::array<V3D,4>& posIn){
    posOut = (timeIn+1.0)*posIn[2];
  }
  void jacTransform(MXD& J,const double& timeIn,const std::array<V3D,4>& posIn){
//    J.resize(n_,m_);
//    J.block<3,1>(0,0) = b
//
//
//
//        getjac<0,0>(MXD* J) =


//    Eigen::Matrix<double,5,1> v;
//    v << 1,2,3,4,5;
//    out.setZero();
//    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()).setIdentity();
//    out.block(posOut_->getIndex(),posIn_->getIndex(),posOut_->getDim(),posIn_->getDim()) *= (timeOffset_->get(in)+2.0);
//    out.block(posOut_->getIndex(),timeOffset_->getIndex(),posOut_->getDim(),timeOffset_->getDim()) = posIn_->get(in)+v;
  }
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
  TransformationExample t;
  State* s1a = t.inputDefinition()->newState();
  State* s1b = t.inputDefinition()->newState();
  s1a->init();
  s1a->print();

  // Boxplus and boxminus
  Eigen::VectorXd v(s1a->getDim());
  v.setZero();
  for(int i = 0; i < s1a->getDim(); i++){
    v(i) = i;
  }
  s1a->boxplus(v,s1b);
  s1b->print();
  s1a->boxminus(s1b,v);
  std::cout << v.transpose() << std::endl;
//
//  // Transformation
//  State* s2 = t.outputDefinition()->newState();
//  MXD J(t.outputDefinition()->getDim(),t.inputDefinition()->getDim());
//  t.jac(s1a,J);
//  std::cout << J << std::endl;
//  MXD JFD(t.outputDefinition()->getDim(),t.inputDefinition()->getDim());
//  t.jacFD(s1a,JFD,1e-8);
//  std::cout << JFD << std::endl;
//
//  delete s1a;
//  delete s1b;
////  delete s2;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
