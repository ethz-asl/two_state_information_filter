#include "gtest/gtest.h"

#include "../include/generalized_information_filter/element-vector.h"
#include "generalized_information_filter/binary-residual.h"
#include "generalized_information_filter/common.h"
#include "generalized_information_filter/filter.h"
#include "generalized_information_filter/prediction.h"
#include "generalized_information_filter/residuals/imu-prediction.h"
#include "generalized_information_filter/residuals/pose-update.h"
#include "generalized_information_filter/transformation.h"
#include "generalized_information_filter/unary-update.h"

using namespace GIF;

class TransformationExample : public Transformation<ElementPack<V3D>,
    ElementPack<double, std::array<V3D, 4>>> {
 public:
  TransformationExample()
      : mtTransformation({"pos"}, {"tim", "sta"}) {
  }

  virtual ~TransformationExample() {}

  void evalTransform(V3D& posOut, const double& timeIn,
                     const std::array<V3D, 4>& posIn) const {
    posOut = (timeIn + 1.0) * (posIn[2] + V3D(1, 2, 3));
  }
  void jacTransform(MXD& J, const double& timeIn,
                    const std::array<V3D, 4>& posIn) const {
    J.setZero();
    setJacBlock<0, 0>(J, V3D(1, 2, 3));
    Eigen::Matrix<double, 3, 12> J2;
    J2.setZero();
    J2.block<3, 3>(0, 6) = M3D::Identity();
    setJacBlock<0, 1>(J, J2);
  }
};

class EmptyMeas : public ElementVector {
 public:
  EmptyMeas(): ElementVector(std::shared_ptr<ElementVectorDefinition>(new ElementVectorDefinition())){};
};

class BinaryRedidualVelocity : public BinaryResidual<ElementPack<V3D>,
    ElementPack<V3D, V3D>, ElementPack<V3D>, ElementPack<V3D>, EmptyMeas> {
 public:
  BinaryRedidualVelocity()
      : mtBinaryRedidual( { "pos" }, { "pos", "vel" }, { "pos" }, { "pos" },
                         false, false, false) {
    dt_ = 0.1;
  }

  virtual ~BinaryRedidualVelocity() {}

  void evalResidualImpl(V3D& posRes, const V3D& posPre, const V3D& velPre,
                        const V3D& posPos, const V3D& posNoi) const {
    posRes = posPre + dt_ * velPre - posPos + posNoi;
  }

  void jacPreImpl(MXD& J, const V3D& posPre, const V3D& velPre,
                  const V3D& posPos, const V3D& posNoi) const {
    J.setZero();
    setJacBlockPre<0, 0>(J, M3D::Identity());
    setJacBlockPre<0, 1>(J, dt_ * M3D::Identity());
  }

  void jacPosImpl(MXD& J, const V3D& posPre, const V3D& velPre,
                  const V3D& posPos, const V3D& posNoi) const {
    J.setZero();
    setJacBlockPos<0, 0>(J, -M3D::Identity());
  }

  void jacNoiImpl(MXD& J, const V3D& posPre, const V3D& velPre,
                  const V3D& posPos, const V3D& posNoi) const {
    J.setZero();
    setJacBlockNoi<0, 0>(J, M3D::Identity());
  }

 protected:
  double dt_;
};

class AccelerometerMeas : public ElementVector {
 public:
  AccelerometerMeas(const V3D& acc = V3D(0, 0, 0))
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(new ElementPack<V3D>( { "acc" }))),
        acc_(ElementVector::GetValue<V3D>("acc")) {
    acc_ = acc;
  }
  V3D& acc_;
};

class BinaryRedidualAccelerometer : public BinaryResidual<ElementPack<V3D>,
    ElementPack<V3D>, ElementPack<V3D>, ElementPack<V3D>, AccelerometerMeas> {
 public:
  BinaryRedidualAccelerometer()
      : mtBinaryRedidual( { "vel" }, { "vel" }, { "vel" }, { "vel" }, false,
                         true, true) {
    dt_ = 0.1;
  }

  virtual ~BinaryRedidualAccelerometer() {
  }

  void evalResidualImpl(V3D& velRes, const V3D& velPre, const V3D& velPos,
                        const V3D& velNoi) const {
    velRes = velPre + dt_ * meas_->acc_ - velPos + velNoi;
  }
  void jacPreImpl(MXD& J, const V3D& velPre, const V3D& velPos,
                  const V3D& velNoi) const {
    J.setZero();
    setJacBlockPre<0, 0>(J, M3D::Identity());
  }
  void jacPosImpl(MXD& J, const V3D& velPre, const V3D& velPos,
                  const V3D& velNoi) const {
    J.setZero();
    setJacBlockPos<0, 0>(J, -M3D::Identity());
  }
  void jacNoiImpl(MXD& J, const V3D& velPre, const V3D& velPos,
                  const V3D& velNoi) const {
    J.setZero();
    setJacBlockNoi<0, 0>(J, M3D::Identity());
  }

 protected:
  double dt_;
};

class PredictionAccelerometer : public Prediction<ElementPack<V3D>,
    ElementPack<V3D>, AccelerometerMeas> {
 public:
  PredictionAccelerometer()
      : mtPrediction( { "vel" }, { "vel" }) {
    dt_ = 0.1;
  }

  virtual ~PredictionAccelerometer() {}

  void evalPredictionImpl(V3D& velPos, const V3D& velPre,
                          const V3D& velNoi) const {
    velPos = velPre + dt_ * meas_->acc_ + velNoi;
  }
  void jacPrePredictionImpl(MXD& J, const V3D& velPre,
                            const V3D& velNoi) const {
    J.setZero();
    setJacBlockPre<0, 0>(J, M3D::Identity());
  }
  void jacNoiPredictionImpl(MXD& J, const V3D& velPre,
                            const V3D& velNoi) const {
    J.setZero();
    setJacBlockNoi<0, 0>(J, M3D::Identity());
  }

 protected:
  double dt_;
};

// The fixture for testing class ScalarState
class NewStateTest : public virtual ::testing::Test {
 protected:
  NewStateTest()
      : covMat_(1, 1) {}
  virtual ~NewStateTest() {}
  MXD covMat_;
};

// Test constructors
TEST_F(NewStateTest, constructor) {
  TransformationExample t;
  std::shared_ptr<ElementVector> s1a(new ElementVector(t.inputDefinition()));
  std::shared_ptr<ElementVector> s1b(new ElementVector(t.inputDefinition()));
  s1a->SetIdentity();
  s1a->Print();

  // Boxplus and BoxMinus
  Eigen::VectorXd v(s1a->GetDimension());
  v.setZero();
  for (int i = 0; i < s1a->GetDimension(); i++) {
    v(i) = i;
  }
  s1a->BoxPlus(v, s1b);
  s1b->Print();
  s1a->BoxMinus(s1b, v);
  std::cout << v.transpose() << std::endl;

  // Jacobian
  MXD J;
  t.jacFD(J, s1a);
  std::cout << J << std::endl;

  // Transformation
  std::shared_ptr<ElementVector> s2(new ElementVector(t.outputDefinition()));
  MXD P1(s1a->GetDimension(), s1a->GetDimension());
  MXD P2(s2->GetDimension(), s2->GetDimension());
  t.transformState(s2, s1a);
  t.transformCovMat(P2, s1a, P1);
  t.testJac(s1a);

  // Velocity Residual
  std::shared_ptr<BinaryRedidualVelocity> velRes(new BinaryRedidualVelocity());
  std::shared_ptr<ElementVector> pre(new ElementVector(velRes->preDefinition()));
  pre->SetIdentity();
  std::shared_ptr<ElementVector> pos(new ElementVector(velRes->posDefinition()));
  pos->SetIdentity();
  std::shared_ptr<ElementVector> noi(new ElementVector(velRes->noiDefinition()));
  noi->SetIdentity();
  velRes->testJacs(pre, pos, noi);

  // Accelerometer Residual
  std::shared_ptr<BinaryRedidualAccelerometer> accRes(
      new BinaryRedidualAccelerometer());

  // Filter
  Filter f;
  f.addRes(velRes);
  f.addRes(accRes);
  std::shared_ptr<ElementVector> preState(new ElementVector(f.stateDefinition()));
  preState->SetIdentity();
  preState->GetValue < V3D > ("pos") = V3D(1, 2, 3);
  preState->Print();
  std::shared_ptr<ElementVector> posState(new ElementVector(f.stateDefinition()));
  posState->SetIdentity();
  posState->Print();
  f.evalRes(preState, posState);


  // Test measurements
  std::shared_ptr<EmptyMeas> eptMeas(new EmptyMeas);
  TimePoint start = Clock::now();
  f.init(start+fromSec(0.00));
  f.addMeas(0, eptMeas,start+fromSec(-0.1));
  f.addMeas(0, eptMeas,start+fromSec(0.0));
  f.addMeas(0, eptMeas,start+fromSec(0.2));
  f.addMeas(0, eptMeas,start+fromSec(0.3));
  f.addMeas(0, eptMeas,start+fromSec(0.4));
  f.addMeas(1, std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(-0.1,0.0,0.0))),start+fromSec(-0.1));
  f.addMeas(1, std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.0,0.0,0.0))),start+fromSec(0.0));
  f.addMeas(1, std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.1,0.0,0.0))),start+fromSec(0.1));
  f.addMeas(1, std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.4,0.0,0.0))),start+fromSec(0.3));
  f.addMeas(1, std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.3,0.0,0.0))),start+fromSec(0.5));

  f.update();

  f.update();

  // Prediction Accelerometer
  std::shared_ptr<PredictionAccelerometer> accPre(
      new PredictionAccelerometer());
  std::shared_ptr<ElementVector> preAcc(new ElementVector(accPre->preDefinition()));
  std::shared_ptr<ElementVector> posAcc(new ElementVector(accPre->preDefinition()));
  std::shared_ptr<ElementVector> noiAcc(new ElementVector(accPre->noiDefinition()));
  preAcc->SetIdentity();
  posAcc->SetIdentity();
  noiAcc->SetIdentity();
  accPre->testJacs(preAcc, posAcc, noiAcc);

  // Test measurements
  Filter f2;
  f2.addRes(velRes);
  f2.addRes(accPre);
  f2.init(start+fromSec(0.00));
  f2.addMeas(0,eptMeas,start+fromSec(-0.1));
  f2.addMeas(0,eptMeas,start+fromSec(0.0));
  f2.addMeas(0,eptMeas,start+fromSec(0.2));
  f2.addMeas(0,eptMeas,start+fromSec(0.3));
  f2.addMeas(0,eptMeas,start+fromSec(0.4));
  f2.addMeas(1,std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(-0.1,0.0,0.0))),start+fromSec(-0.1));
  f2.addMeas(1,std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.0,0.0,0.0))),start+fromSec(0.0));
  f2.addMeas(1,std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.1,0.0,0.0))),start+fromSec(0.1));
  f2.addMeas(1,std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.4,0.0,0.0))),start+fromSec(0.3));
  f2.addMeas(1,std::shared_ptr<AccelerometerMeas>(
      new AccelerometerMeas(V3D(0.3,0.0,0.0))),start+fromSec(0.5));
  f2.update();
  f2.update();

  // Test IMU + Pose filter
  int s = 0;
  std::shared_ptr<IMUPrediction> imuPre(new IMUPrediction());
  imuPre->getR() = 1e-8 * imuPre->getR();
  imuPre->testJacs(s);
  std::shared_ptr<PoseUpdate> poseUpd(new PoseUpdate());
  poseUpd->getR() = 1e-8 * poseUpd->getR();
  poseUpd->testJacs(s);

  Filter imuPoseFilter;
  int imuPreInd = imuPoseFilter.addRes(imuPre);
  int poseUpdInd = imuPoseFilter.addRes(poseUpd);
  imuPoseFilter.init(start);
  imuPoseFilter.addMeas(imuPreInd, std::shared_ptr<IMUMeas>(
      new IMUMeas(V3D(0.0, 0.0, 0.0), V3D(0.0, 0.0, 9.81))), start);
  for (int i = 1; i <= 10; i++) {
    imuPoseFilter.addMeas(imuPreInd, std::shared_ptr<IMUMeas>(
        new IMUMeas(V3D(0.3, 0.0, 0.1), V3D(0.0, 0.2, 9.81))),
        start + fromSec(i));
    imuPoseFilter.addMeas(poseUpdInd, std::shared_ptr<PoseMeas>(
        new PoseMeas(V3D(0.0, 0.0, 0.0), QPD(1.0, 0.0, 0.0, 0.0))),
        start + fromSec(i));
  }
  imuPoseFilter.update();

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
