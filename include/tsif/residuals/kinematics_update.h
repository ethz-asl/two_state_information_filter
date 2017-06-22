#pragma once

//base class
#include "tsif/residual.h"

namespace tsif {

//measurement class to hold base to foot vector and contact flag
class ContactPointMeasurement : public ElementVector<Element<Vec3,0>, Element<bool,1>>{
 public:
  //constructors with and without arguments
  ContactPointMeasurement(): ElementVector<Element<Vec3,0>,Element<bool,1>>(Vec3(.0,.0,.0),false) {}
  ContactPointMeasurement(const Vec3& B_Ns,const bool& Nlambda):ElementVector<Element<Vec3,0>, Element<bool,1>>(B_Ns,Nlambda){}
  //functions to get contact point
  const Vec3& GetContactPoint() const {return Get<0>();}
  Vec3& GetContactPoint() {return Get<0>();}
  //function to get contact flag
  const bool& GetContactFlag() const {return Get<1>();}
  bool& GetContactFlag() {return Get<1>();}
};

//shortcut to make the templating more convenient
template<int Y,//contact point innovation
         int I_R_IB,//position base to odom in odom
         int PHI_IB,//orientation odom to base
         int I_P>//right front contact point in odom
using KinematicsUpdateBase = Residual<ElementVector<Element<Vec3,Y>>,//innovation
                                      ElementVector<>,//previous state is empty
                                      ElementVector<Element<Vec3,I_R_IB>,Element<Quat,PHI_IB>,Element<Vec3,I_P>>,//current state
                                      ContactPointMeasurement>;//measurement

//kinematics update residual comparing measured contact points calculated from joint state measurements to contact points in the filter state
template<int I_R_IB,int PHI_IB,int I_P>
class KinematicsUpdate : public KinematicsUpdateBase<0,I_R_IB,PHI_IB,I_P> {
 public:
  //typedef and types from base for convenience
  //base type
  typedef KinematicsUpdateBase<0,I_R_IB,PHI_IB,I_P> Base;
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
  //constructor setting characteristics of the residual and default weights
  KinematicsUpdate(): Base (false,//splittable
                            false,//mergeable
                            false),//mandatory
                            w_contact_(1.), w_swing_(1.) {};
  //function evaluating the innovations
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
   //compare measured contact points to filter state
		 out. template Get<0>() = meas_->GetContactPoint()-cur.template Get<PHI_IB>().inverse().toRotationMatrix()*(cur.template Get<I_P>()-cur.template Get<I_R_IB>());
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
   //derivative of contact point innovation wrt position
   this->template SetJacCur<0,I_R_IB>(J,cur,C_BI);
   //derivative of contact point innovation wrt orientation
   this->template SetJacCur<0,PHI_IB>(J,cur,-SSM(C_BI*(cur.template Get<I_P>()-cur.template Get<I_R_IB>()))*C_BI);
   //derivative of contact point innovation wrt contact point
   this->template SetJacCur<0,I_P>(J,cur,-C_BI);
   return 0;    
  }
  //function returning the weighting factor used in the TSIF optimization
  double GetWeight(){
    if (meas_->GetContactFlag()) return w_contact_/sqrt(dt_);
    else return w_swing_/sqrt(dt_);
  }
  //weights for contact points during contact and swing
  double w_contact_;
  double w_swing_;
};

} // namespace tsif

