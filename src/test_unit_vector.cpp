#include "tsif/unit_vector.h"

using namespace tsif;

int main(int argc, char** argv){
  tsif::UnitVector uv(Vec3(1,2,3));
  std::cout << uv.GetVec().transpose() << std::endl;
  std::cout << Vec3(1,2,3).normalized().transpose() << std::endl;
  uv.SetIdentity();
  std::cout << uv.GetVec().transpose() << std::endl;
  uv.SetFromVector(Vec3(1,2,3)*1e-12);
  std::cout << uv.GetVec().transpose() << std::endl;
  uv.SetFromVector(Vec3(0,0,0));
  std::cout << uv.GetVec().transpose() << std::endl;
  uv.SetFromVector(Vec3(0,0,-1));
  std::cout << uv.GetVec().transpose() << std::endl;

  tsif::UnitVector uv0(Vec3(0.1,0,1));
  tsif::UnitVector uv1(Vec3(0.2,0,1));
  std::cout << uv1.GetVec().transpose() << std::endl;
  Vec<2> unitVectorDif;
  uv1.Boxminus(uv0,unitVectorDif);
  std::cout << unitVectorDif.transpose() << std::endl;
  uv1.SetIdentity();
  std::cout << uv1.GetVec().transpose() << std::endl;
  uv0.Boxplus(unitVectorDif,uv1);
  std::cout << uv1.GetVec().transpose() << std::endl;

  tsif::UnitVector test_uv0, test_uv0_pert, test_uv1, test_uv1_pert;
  test_uv0.SetRandom();
  test_uv1.SetRandom();
  Vec<2> dir = tsif::NormalRandomNumberGenerator::Instance().GetVec<2>();
  Vec<2> ref,pert;
  test_uv1.Boxminus(test_uv0,ref);
  test_uv0.Boxplus(dir*1e-8,test_uv0_pert);
  test_uv1.Boxplus(dir*1e-8,test_uv1_pert);
  test_uv1.Boxminus(test_uv0_pert,pert);
  std::cout << "BoxminusJacRef" << std::endl;
  std::cout << (pert-ref).transpose()*1e8 << std::endl;
  Mat<2> J;
  test_uv1.BoxminusJacRef(test_uv0,J);
  std::cout << (J*dir).transpose() << std::endl;
  test_uv1_pert.Boxminus(test_uv0,pert);
  std::cout << "BoxminusJacInp" << std::endl;
  std::cout << (pert-ref).transpose()*1e8 << std::endl;
  test_uv1.BoxminusJacInp(test_uv0,J);
  std::cout << (J*dir).transpose() << std::endl;

  ref = tsif::NormalRandomNumberGenerator::Instance().GetVec<2>();
  test_uv0.SetRandom();
  test_uv0.Boxplus(ref,test_uv1);
  test_uv0.Boxplus(ref+dir*1e-8,test_uv1_pert);
  test_uv1_pert.Boxminus(test_uv1,pert);
  const Vec<2> der1 = pert*1e8;
  test_uv0.BoxplusJacVec(ref,J);
  const Vec<2> der2 = J*dir;
  std::cout << "BoxplusJavVec" << std::endl;
  std::cout << der1.transpose() << std::endl;
  std::cout << der2.transpose() << std::endl;
  test_uv0.Boxplus(ref,test_uv1);
  test_uv0.Boxplus(dir*1e-8,test_uv0_pert);
  test_uv0_pert.Boxplus(ref,test_uv1_pert);
  test_uv1_pert.Boxminus(test_uv1,pert);
  test_uv0.BoxplusJacInp(ref,J);
  std::cout << "BoxplusJacInp" << std::endl;
  std::cout << (pert).transpose()*1e8 << std::endl;
  std::cout << (J*dir).transpose() << std::endl;

  return 0;
}
