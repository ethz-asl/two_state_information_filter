#include "tsif/element_vector.h"

using namespace tsif;

int main(int argc, char** argv){

  typedef ElementVector<Element<Vec3,0>,Element<double,1>> ElementVectorA;
  typedef ElementVector<Element<double,2>,Element<double,1>> ElementVectorB;
  typedef ElementVector<Element<Quat,3>,Element<std::array<Vec<2>,2>,4>> ElementVectorC;
  typedef typename MergeTrait<ElementVectorA,ElementVectorB,ElementVectorC>::Type ElementVectorD;

  ElementVectorD vec;
  vec.Get<0>() = Vec3(1,-2,3);
  vec.Get<1>() = 1.5;
  std::cout << "==== ElementVector ====\n" << vec.Print() << "=======================" << std::endl;
  vec.SetIdentity();
  std::cout << "==== ElementVector ====\n" << vec.Print() << "=======================" << std::endl;
  vec.SetRandom();
  std::cout << "==== ElementVector ====\n" << vec.Print() << "=======================" << std::endl;

  std::cout << ElementVectorD::Start(0) << std::endl;
  std::cout << ElementVectorD::Start(1) << std::endl;
  std::cout << ElementVectorD::Start(2) << std::endl;
  std::cout << ElementVectorD::Start(3) << std::endl;
  std::cout << ElementVectorD::Start(4) << std::endl;
  std::cout << ElementVectorD::Dim() << std::endl;

  ElementVectorD vec2;
  vec2.SetRandom();
  std::cout << "==== ElementVector ====\n" << vec2.Print() << "=======================" << std::endl;
  Vec<ElementVectorD::Dim()> dif;
  vec2.Boxminus(vec,dif);
  std::cout << dif.transpose() << std::endl;
  vec2.SetIdentity();
  std::cout << "==== ElementVector ====\n" << vec2.Print() << "=======================" << std::endl;
  vec.Boxplus(dif,vec2);
  std::cout << "==== ElementVector ====\n" << vec2.Print() << "=======================" << std::endl;

  ElementVector<Element<Vec3,0>> cur;
  cur.SetRandom();
  ElementVector<Element<Vec3,0>,Element<Vec3,1>> pre;
  pre.SetRandom();
  std::cout << pre.Print() << std::endl;
  ElementVectorRef<Element<Vec3,1>> pre_ref(pre);
  std::cout << pre_ref.Print() << std::endl;

  return 0;
}
