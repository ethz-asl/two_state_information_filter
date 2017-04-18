#ifndef TSIF_MODEL_H_
#define TSIF_MODEL_H_

#include "tsif/utils/common.h"
#include "element_vector.h"

namespace tsif{

template<typename Derived, typename Out, typename... Ins>
class Model{
 public:
  ~Model(){};
  int Eval(typename Out::Ref out, const std::tuple<typename Ins::CRef...> ins){
    return static_cast<Derived&>(*this).template EvalImpl(out,ins);
  }
  template<int N>
  int Jac(MatRefX J, const std::tuple<typename Ins::CRef...> ins){
    return static_cast<Derived&>(*this).template JacImpl<N>(J,ins);
  }
  template<int N, int C>
  int JacFindif(double d, MatRefX J, const std::tuple<typename Ins::CRef...> insRef){
    const int outDim = Out::Dim();
    const int I = std::tuple_element<N,std::tuple<Ins...>>::type::template GetId<C>();
    const int inDim = std::tuple_element<N,std::tuple<Ins...>>::type::template GetElementDim<I>();
    Out outRef,outDis;
    Eval(outRef,insRef);
    std::tuple<Ins...> insDis = insRef;
    Vec<inDim> inDif;
    Vec<outDim> outDif;
    for(unsigned int j=0;j<inDim;j++){
      inDif.setZero();
      inDif(j) = d;
      std::get<N>(insRef).template GetElement<I>().Boxplus(inDif,std::get<N>(insDis).template GetElement<I>());
      Eval(outDis,insDis);
      outDis.Boxminus(outRef,outDif);
      J.col(std::get<N>(insRef).Start(I)+j) = outDif/d;
    }
    return 0;
  }
  template<int N, int C = 0, typename std::enable_if<(C < std::tuple_element<N,std::tuple<Ins...>>::type::kN)>::type* = nullptr>
  int JacFindifFull(double d, MatRefX J, const std::tuple<typename Ins::CRef...> insRef){
    JacFindif<N,C>(d,J,insRef);
    JacFindifFull<N,C+1>(d,J,insRef);
    return 0;
  }
  template<int N, int C = 0, typename std::enable_if<(C >= std::tuple_element<N,std::tuple<Ins...>>::type::kN)>::type* = nullptr>
  int JacFindifFull(double d, MatRefX J, const std::tuple<typename Ins::CRef...> insRef){
    return 0;
  }
  template<int N>
  int JacTest(double th, double d, const std::tuple<typename Ins::CRef...> ins){
    const int outDim = Out::Dim();
    const int inDim = std::get<N>(ins).Dim();
    if(outDim > 0 & inDim > 0){
      MatX J(outDim,inDim);
      J.setZero();
      MatX J_FD(outDim,inDim);
      J_FD.setZero();
      Jac<N>(J,ins);
      JacFindifFull<N>(d,J_FD,ins);
      TSIF_LOG("Analytical Jacobian:\n" << J);
      TSIF_LOG("Numerical Jacobian:\n" << J_FD);
      typename MatX::Index maxRow, maxCol = 0;
      const double r = (J-J_FD).array().abs().maxCoeff(&maxRow, &maxCol);
      if(r>th){
        std::cout << "\033[31m==== Model jacInput (" << N << ") Test failed: " << r
                  << " is larger than " << th << " at row "
                  << maxRow << " and col " << maxCol << " ===="
                  << "  " << J(maxRow,maxCol) << " vs " << J_FD(maxRow,maxCol)
                  << "\033[0m" << std::endl;
        return 1;
      } else{
        std::cout << "\033[32m==== Test successful (" << r << ") ====\033[0m" << std::endl;
        return 0;
      }
    } else {
      std::cout << "\033[32m==== Test successful ( dimension 0 ) ====\033[0m" << std::endl;
    }
  }
};

} // namespace tsif

#endif  // TSIF_MODEL_H_
