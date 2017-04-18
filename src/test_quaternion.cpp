#include "tsif/utils/common.h"
#include "tsif/element.h"

using namespace tsif;

int main(int argc, char** argv){
  Vec3 vec = NormalRandomNumberGenerator::Instance().GetVec<3>();
  Quat q(1,0,0,0);
  ElementTraits<Quat>::SetRandom(q);
  Quat p(1,0,0,0);
  ElementTraits<Quat>::SetRandom(p);
  std::cout << (q.toRotationMatrix()*p.toRotationMatrix()*vec).transpose() << std::endl;
  std::cout << ((q*p).toRotationMatrix()*vec).transpose() << std::endl;
  std::cout << q.coeffs().transpose() << std::endl;
  std::cout << Exp(Log(q)).coeffs().transpose() << std::endl;
  std::cout << RotMat(Log(q)) << std::endl;
  std::cout << q.toRotationMatrix() << std::endl;

  for(int i=0;i<10000;i++){
    Quat q1;
    ElementTraits<Quat>::SetRandom(q1);
    Quat q2 = Exp(Log(q1));
    if(Boxminus(q1,q2).norm() > 1e-8){
      std::cout << "Fault" << std::endl;
    }
  }

  return 0;
}
