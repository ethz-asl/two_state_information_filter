#include "gtest/gtest.h"

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"
#include "generalized_information_filter/Transformation.hpp"
#include "generalized_information_filter/BinaryResidual.hpp"
#include "generalized_information_filter/Filter.hpp"

using namespace GIF;

class TransformationExample: public Transformation<ElementPack<V3D>,ElementPack<double, std::array<V3D,4>>>{
 public:
  TransformationExample(): mtTransformation({"pos"},{"tim","sta"}){};
  virtual ~TransformationExample(){};
  void evalTransform(V3D& posOut,const double& timeIn,const std::array<V3D,4>& posIn) const{
    posOut = (timeIn+1.0)*(posIn[2]+V3D(1,2,3));
  }
  void jacTransform(MXD& J,const double& timeIn,const std::array<V3D,4>& posIn) const{
    J.setZero();
    setJacBlock<0,0>(J,V3D(1,2,3));
    Eigen::Matrix<double,3,12> J2;
    J2.setZero();
    J2.block<3,3>(0,6) = M3D::Identity();
    setJacBlock<0,1>(J,J2);
  }
};

class AccelerometerMeas: public BinaryMeasurementBase{
 public:
  AccelerometerMeas(): BinaryMeasurementBase(ElementPack<V3D>::makeStateDefinition({"acc"})),
      acc_(State::getValue<V3D>("acc")){};
  V3D& acc_;
};

class EmptyMeas: public BinaryMeasurementBase{
 public:
  EmptyMeas(): BinaryMeasurementBase(std::shared_ptr<StateDefinition>(new StateDefinition())){};
};

class BinaryRedidualVelocity: public BinaryResidual<ElementPack<V3D>,ElementPack<V3D,V3D>,ElementPack<V3D>,ElementPack<V3D>,EmptyMeas>{
 public:
  BinaryRedidualVelocity(): mtBinaryRedidual({"pos"},{"pos","vel"},{"pos"},{"pos"}){
    dt_ = 0.1;
  };
  virtual ~BinaryRedidualVelocity(){};
  void evalResidual(V3D& posRes,const V3D& posPre,const V3D& velPre,const V3D& posPos,const V3D& posNoi) const{
    posRes = posPre + dt_*velPre - posPos + posNoi;
  }
  void jacPre(MXD& J,const V3D& posPre,const V3D& velPre,const V3D& posPos,const V3D& posNoi) const{
    J.setZero();
    setJacBlockPre<0,0>(J,M3D::Identity());
    setJacBlockPre<0,1>(J,dt_*M3D::Identity());
  }
  void jacPos(MXD& J,const V3D& posPre,const V3D& velPre,const V3D& posPos,const V3D& posNoi) const{
    J.setZero();
    setJacBlockPre<0,0>(J,-M3D::Identity());
  }
  void jacNoi(MXD& J,const V3D& posPre,const V3D& velPre,const V3D& posPos,const V3D& posNoi) const{
    J.setZero();
    setJacBlockPre<0,0>(J,M3D::Identity());
  }

 protected:
  double dt_;
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
  std::shared_ptr<State> s1a(new State(t.inputDefinition()));
  std::shared_ptr<State> s1b(new State(t.inputDefinition()));
  s1a->setIdentity();
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

  // Jacobian
  MXD J;
  t.jacFD(J,s1a);
  std::cout << J << std::endl;

  // Transformation
  std::shared_ptr<State> s2(new State(t.outputDefinition()));
  MXD P1(s1a->getDim(),s1a->getDim());
  MXD P2(s2->getDim(),s2->getDim());
  t.transformState(s2,s1a);
  t.transformCovMat(P2,s1a,P1);
  t.testJac(s1a);

  // Residual
  std::shared_ptr<BinaryRedidualVelocity> velRes(new BinaryRedidualVelocity());
  std::shared_ptr<State> pre(new State(velRes->preDefinition()));
  pre->setIdentity();
  std::shared_ptr<State> pos(new State(velRes->posDefinition()));
  pos->setIdentity();
  std::shared_ptr<State> noi(new State(velRes->noiDefinition()));
  noi->setIdentity();
  velRes->testJacs(pre,pos,noi);

  // Filter
  Filter f;
  f.addRes(velRes);
  std::shared_ptr<State> preState(new State(f.stateDefinition()));
  preState->setIdentity();
  preState->getValue<V3D>("pos") = V3D(1,2,3);
  preState->print();
  std::shared_ptr<State> posState(new State(f.stateDefinition()));
  posState->setIdentity();
  posState->print();
  f.evalRes(preState,posState);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
