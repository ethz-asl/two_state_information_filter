#ifndef TSIF_UNIT_VECTOR_H_
#define TSIF_UNIT_VECTOR_H_

#include "tsif/utils/common.h"

namespace tsif{

class UnitVector{
 private:
  Quat q_;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  UnitVector(){}
  UnitVector(const UnitVector& other): q_(other.GetQuat()){}
  UnitVector(const Quat& q): q_(q){}
  UnitVector(const Vec3& v){
    SetFromVector(v);
  }
  const Quat& GetQuat() const{
    return q_;
  }
  void SetQuat(const Quat& q){
    q_ = q;
  }
  Vec3 GetVec() const{
    return q_.toRotationMatrix()*Vec3(0,0,1);
  }
  Vec3 GetPerp1() const{
    return q_.toRotationMatrix()*Vec3(1,0,0);
  }
  Vec3 GetPerp2() const{
    return q_.toRotationMatrix()*Vec3(0,1,0);
  }
  void SetFromVector(const Vec3& vec){
    q_ = vec.norm() < 1e-12 ? Quat::Identity() : Quat::FromTwoVectors(Vec3(0,0,1),vec);
  }
  void SetIdentity(){
    q_.setIdentity();
  }
  void SetRandom(){
    SetFromVector(NormalRandomNumberGenerator::Instance().GetVec<3>());
  }
  Eigen::Matrix<double,3,2> GetM() const {
    Eigen::Matrix<double,3,2> M;
    M.col(0) = -GetPerp2();
    M.col(1) = GetPerp1();
    return M;
  }
  Eigen::Matrix<double,3,2> GetN() const {
    Eigen::Matrix<double,3,2> M;
    M.col(0) = GetPerp1();
    M.col(1) = GetPerp2();
    return M;
  }
  void Boxplus(const VecCRef<2>& dif, UnitVector& out) const{
    out.SetQuat(Exp(dif(0)*GetPerp1()+dif(1)*GetPerp2())*q_);
  }
  void Boxminus(const UnitVector& ref, VecRef<2> dif) const{
    dif = ref.GetN().transpose()*Log(Quat::FromTwoVectors(ref.GetVec(),GetVec()));
  }
  void BoxminusJacRef(const UnitVector& ref, MatRef<2> J) const{
    J = ref.GetN().transpose()*FromTwoVectorsJac(ref.GetVec(),GetVec())*ref.GetM();
  }
  void BoxminusJacInp(const UnitVector& ref, MatRef<2> J) const{
    J = -ref.GetN().transpose()*FromTwoVectorsJac(GetVec(),ref.GetVec())*GetM();
  }
  void BoxplusJacVec(const VecCRef<2>& dif, MatRef<2> J) const{
    UnitVector out;
    Boxplus(dif,out);
    J = out.GetN().transpose()*GammaMat(dif(0)*GetPerp1()+dif(1)*GetPerp2())*GetN();
  }
  void BoxplusJacInp(const VecCRef<2>& dif, MatRef<2> J) const{
    UnitVector out;
    Boxplus(dif,out);
    J = out.GetN().transpose()*GetN();
  }
};

} // namespace tsif

#endif  // TSIF_UNIT_VECTOR_H_
