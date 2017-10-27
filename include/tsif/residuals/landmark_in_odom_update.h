#pragma once

//base class
#include "tsif/residual.h"
#include "tsif/utils/common.h"

namespace tsif {

//measurement class to hold base to landmark vector in odom frame
class LandmarkInOdomMeasurement : public ElementVector<Element<Vec3,0>>{
 public:
  //constructors with and without arguments
  LandmarkInOdomMeasurement(): ElementVector<Element<Vec3,0>>(Vec3(.0,.0,.0)) {}
  LandmarkInOdomMeasurement(const Vec3& B_Ns):ElementVector<Element<Vec3,0>>(B_Ns){}
  //functions to measurement
  const Vec3& GetMeasurement() const {return Get<0>();}
  Vec3& GetMeasurement() {return Get<0>();}
};

//shortcut to make the templating more convenient
template<int Y,//contact point innovation
         int I_R_IB,//position base to odom in odom
         int PHI_IB,//orientation odom to base
         int I_P>//right front contact point in odom
using LandmarkInOdomUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,//innovation
                                          ElementVector<>,//previous state is empty
                                          ElementVector<Element<Vec3,I_R_IB>,Element<Quat,PHI_IB>,Element<Vec3,I_P>>,//current state
                                          LandmarkInOdomMeasurement>;//measurement

//landmark update residual comparing measured contact points calculated from joint state measurements to contact points in the filter state
template<int I_R_IB,int PHI_IB,int I_P>
class LandmarkInOdomUpdate : public LandmarkInOdomUpdateBase<0,I_R_IB,PHI_IB,I_P> {
 public:
  //typedef and types from base for convenience
  //base type
  typedef LandmarkInOdomUpdateBase<0,I_R_IB,PHI_IB,I_P> Base;
  //innovation vector type
  using typename Base::Output;
  //previous state vector type
  using typename Base::Previous;
  //current state vector type
  using typename Base::Current;
  //time increment from base
  using Base::dt_;
  //measurement pointer from base
  using Base::meas_;
  //weight
  using Base::w_;
  //parameter for the huber loss function
  double huber_threshold_;
  //constructor setting characteristics of the residual and default weights
  LandmarkInOdomUpdate(): Base (false,//splittable
                                false,//mergeable
                                false),//mandatory
                          huber_threshold_(1000.){};
  //function evaluating the innovations
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
   //compare measured contact points to filter state
		 out. template Get<0>() = meas_->GetMeasurement()-cur.template Get<PHI_IB>().inverse().toRotationMatrix()*(cur.template Get<I_P>()-cur.template Get<I_R_IB>());
   return 0;
  }
  //function evaluating the jacobian of the the innovation wrt the previous state
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //has dimension 0
   return 0;    
  }
  //function evaluating the jacobian of the the innovation wrt the current state
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
   //set jacobian zero
   J.setZero();
   //rotation matrix from odom to base frame from orientation state
   const Mat3 C_BI=cur.template Get<PHI_IB>().inverse().toRotationMatrix();
   //set nonzero blocks of the jacobian
   //derivative of landmark innovation wrt position
   this->template SetJacCur<0,I_R_IB>(J,cur,C_BI);
   //derivative of landmark innovation wrt orientation
   this->template SetJacCur<0,PHI_IB>(J,cur,-SSM(C_BI*(cur.template Get<I_P>()-cur.template Get<I_R_IB>()))*C_BI);
   //derivative of landmark innovation wrt landmark
   this->template SetJacCur<0,I_P>(J,cur,-C_BI);
   return 0;    
  }
  //function to scale the innovation and jacobians with the appropriate weighting factor
  virtual void AddNoise(typename Output::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    //compute weight using huber loss function
    double weight = w_;
    const double norm = out.template Get<0>().norm();
    //std::cout << "||y_lmk||_2 = " << norm << std::endl;
    if(norm > huber_threshold_) weight *= sqrt(2 * huber_threshold_ * (norm - 0.5 * huber_threshold_)/(norm*norm));
    //scale the innovation and jacobians
    out.Scale(weight);
    J_pre *= weight;
    J_cur *= weight;
  }
};
} // namespace tsif

