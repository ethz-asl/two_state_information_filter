#ifndef GIF_KINEMATICSMODEL_H_
#define GIF_KINEMATICSMODEL_H_

#include <Eigen/Dense>

namespace GIF {

template<int NumLeg, int NumDof>
class LeggedRobotModel{
 public:
  static constexpr double g_ = 9.81;
  static constexpr int kNumLeg = NumLeg;
  static constexpr int kNumDof = NumDof;

  LeggedRobotModel(){};
  virtual ~LeggedRobotModel(){};

  // TODO: cleanup
  virtual Eigen::Vector3d forwardKinematicsBaseToFootInBaseFrame(Eigen::Vector3d angles,unsigned int legId) = 0;
  virtual Eigen::Matrix3d getJacobianTranslationBaseToSegment(Eigen::Vector3d angles,unsigned int legId,unsigned int segmentId) = 0;
  virtual Eigen::Matrix3d getJacobianTranslationBaseToFoot(Eigen::Vector3d angles,unsigned int legId){
    return getJacobianTranslationBaseToSegment(angles,legId,0);
  }
  virtual Eigen::Vector3d forwardKinematicsWorldToFootInWorldFrame(Eigen::Vector3d angles,unsigned int legId) = 0;
  virtual Eigen::Matrix<double,3,18> getJacobianTranslationWorldToSegment(unsigned int legId,unsigned int segmentId) = 0;
  virtual Eigen::Matrix<double,3,18> getJacobianTranslationWorldToFoot(unsigned int legId){
    return getJacobianTranslationWorldToSegment(legId,0);
  }
  virtual void setGeneralizedPositions(const Eigen::VectorXd& genPos, bool update){};
  virtual void setGeneralizedVelocities(const Eigen::VectorXd& genVel, bool update){};
  virtual bool getMassInertiaMatrix(Eigen::MatrixXd& M){return false;};
  virtual bool getNonlinearEffects(Eigen::VectorXd& h){return false;};
  virtual bool getSelectionMatrix(Eigen::MatrixXd& S){return false;};

  // Not used atm
  virtual Eigen::Matrix3d getOrientationBaseToFoot(int legId){return Eigen::Matrix3d::Identity();};
  virtual bool getJacobianFullTranslationBaseToFoot(Eigen::MatrixXd& J, int i){return false;};
};

class LeggedRobotModelExample: public LeggedRobotModel<4,3>{
 public:
  static constexpr double bx = 0.2525;
  static constexpr double by = 0.185;
  static constexpr double lH_0 = -0.0685;
  static constexpr double lT_0 = -0.2;
  static constexpr double lS_0 = -0.235;

  LeggedRobotModelExample(){};

  Eigen::Vector3d forwardKinematicsBaseToFootInBaseFrame(Eigen::Vector3d angles,unsigned int legId){
    double sc0;
    double sc1;
    double sc2;
    double sc3;
    double sc4;
    double sc5;
    sc0 = sin(angles(0));
    sc1 = sin(angles(1));
    sc2 = sin(angles(1)+angles(2));
    sc3 = cos(angles(0));
    sc4 = cos(angles(1));
    sc5 = cos(angles(1)+angles(2));
    Eigen::Vector3d s;
    s(0) = ((legId<2)*2-1)*bx+lT_0*sc1+lS_0*sc2;
    s(1) = -(((int)legId%2)*2-1)*by-sc0*(lH_0+lT_0*sc4+lS_0*sc5);
    s(2) = sc3*(lH_0+lT_0*sc4+lS_0*sc5);
    return s;
  };
	Eigen::Matrix3d getJacobianTranslationBaseToSegment(Eigen::Vector3d angles,unsigned int legId,unsigned int segmentId){
	  double lH;
	  double lT;
	  double lS;

	  switch (segmentId) {
	    case 0:   // foot (contact point to ground)
	      lH = lH_0;
	      lT = lT_0;
	      lS = lS_0;
	      break;
	    case 1:   // CoM of shank-link
	      lH = lH_0;
	      lT = lT_0;
	      lS = lS_0/2;
	      break;
	    case 2:   // CoM of thigh-link
	      lH = lH_0;
	      lT = lT_0/2;
	      lS = 0;
	      break;
	    case 3:   // CoM of hip-link
	      lH = lH_0/2;
	      lT = 0;
	      lS = 0;
	      break;
	    default:
	      std::cout << "no legKinJac assigned to segmentID " << segmentId << "." << std::endl;
	      break;
	  }

	  double sc0;
	  double sc1;
	  double sc2;
	  double sc3;
	  double sc4;
	  double sc5;
	  sc0 = sin(angles(0));
	  sc1 = sin(angles(1));
	  sc2 = sin(angles(1)+angles(2));
	  sc3 = cos(angles(0));
	  sc4 = cos(angles(1));
	  sc5 = cos(angles(1)+angles(2));
	  Eigen::Matrix<double,3,3> J;
	  J.setZero();
	  J(0,1) = lS*sc5+lT*sc4;
	  J(0,2) = lS*sc5;
	  J(1,0) = -sc3*(lH+lT*sc4+lS*sc5);
	  J(1,1) = sc0*(lT*sc1+lS*sc2);
	  J(1,2) = lS*sc0*sc2;
	  J(2,0) = -sc0*(lH+lT*sc4+lS*sc5);
	  J(2,1) = -sc3*(lT*sc1+lS*sc2);
	  J(2,2) = -lS*sc3*sc2;
	  return J;
	};
};

}

#endif /* GIF_KINEMATICSMODEL_H_ */
