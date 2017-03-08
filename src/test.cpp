#include "tsif/residuals/bearing_findif.h"

int main(int argc, char** argv){

  tsif::BearingFindif<0,0,1,2,3,4,5,2> test;
  test.JacPreTest(1e-6,1e-8);
  test.JacCurTest(1e-6,1e-8);


  tsif::UnitVector test_uv0,test_uv1,test_uv0_pert,test_uv1_pert;
  tsif::Vec<2> vec2,dir,pert;
  vec2 = tsif::NormalRandomNumberGenerator::Instance().GetVec<2>();
  dir = tsif::NormalRandomNumberGenerator::Instance().GetVec<2>();
  test_uv0.SetRandom();
  test_uv0.Boxplus(vec2,test_uv1);
  test_uv0.Boxplus(vec2+dir*1e-8,test_uv1_pert);
  test_uv1_pert.Boxminus(test_uv1,pert);
  const tsif::Vec<2> der1 = pert*1e8;
  tsif::Mat<2> J;
  test_uv0.BoxplusJacVec(vec2,J);
  const tsif::Vec<2> der2 = J*dir;
  std::cout << "BoxplusJavVec" << std::endl;
  std::cout << der1.transpose() << std::endl;
  std::cout << der2.transpose() << std::endl;
  test_uv0.Boxplus(vec2,test_uv1);
  test_uv0.Boxplus(dir*1e-8,test_uv0_pert);
  test_uv0_pert.Boxplus(vec2,test_uv1_pert);
  test_uv1_pert.Boxminus(test_uv1,pert);
  test_uv0.BoxplusJacInp(vec2,J);
  std::cout << "BoxplusJacInp" << std::endl;
  std::cout << (pert).transpose()*1e8 << std::endl;
  std::cout << (J*dir).transpose() << std::endl;



  tsif::Vec<2> out2,out2_pert;
  tsif::UnitVector test_uv2;
  tsif::Mat<2> J1,J2;
  const double d = 1e-8;
  // Reference
  test_uv0.Boxplus(vec2,test_uv1);
  test_uv2.Boxminus(test_uv1,out2);
  test_uv0.BoxplusJacInp(vec2,J1);
  test_uv2.BoxminusJacRef(test_uv1,J2);

  // Perturbed
  test_uv0.Boxplus(dir*d,test_uv0_pert);
  test_uv0_pert.Boxplus(vec2,test_uv1_pert);
  test_uv2.Boxminus(test_uv1_pert,out2_pert);

  std::cout << "Jac Test 1" << std::endl;
  tsif::Vec<2> dif;
  test_uv1_pert.Boxminus(test_uv1,dif);
  std::cout << test_uv1_pert.GetQuat().coeffs().transpose() << std::endl;
  std::cout << (dif).transpose()/d << std::endl;
  std::cout << (J1*dir).transpose() << std::endl;

  std::cout << "Jac Test 2" << std::endl;
  std::cout << (out2_pert-out2).transpose()/d << std::endl;
  std::cout << (J2*J1*dir).transpose() << std::endl;



  dir = dif/d;
  test_uv1.Boxplus(dir*d,test_uv1_pert);
  std::cout << test_uv1_pert.GetQuat().coeffs().transpose() << std::endl;
  test_uv2.Boxminus(test_uv1,out2);
  test_uv2.Boxminus(test_uv1_pert,out2_pert);

  std::cout << "Jac Test 3" << std::endl;
  std::cout << (out2_pert-out2).transpose()/d << std::endl;
  std::cout << (J2*dir).transpose() << std::endl;





  return 0;
}
