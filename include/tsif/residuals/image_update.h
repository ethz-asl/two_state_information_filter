#ifndef TSIF_IMAGE_UPDATE_H_
#define TSIF_IMAGE_UPDATE_H_

#include "tsif/utils/common.h"
#include "tsif/residual.h"

namespace tsif{

template<int OUT_RES, int STA_BEA, int N, typename MEAS>
using ImageUpdateBase = Residual<ElementVector<Element<std::array<Vec<2>,N>,OUT_RES>>,
                                 ElementVector<>,
                                 ElementVector<Element<std::array<UnitVector,N>,STA_BEA>>,
                                 MEAS>;

template<int OUT_RES, int STA_BEA, int N, typename MEAS>
class ImageUpdate: public ImageUpdateBase<OUT_RES,STA_BEA,N,MEAS>{ // TODO: rename to bearing update
 public:
  typedef ImageUpdateBase<OUT_RES,STA_BEA,N,MEAS> Base;
  using Base::meas_;
  using Base::w_;
  using Base::GetWeight;
  typedef typename Base::Output Output;
  typedef typename Base::Previous Previous;
  typedef typename Base::Current Current;
  ImageUpdate(): Base(false,false,false){
    for(int i=0;i<N;i++){
      active_[i] = true;
    }
    huberTh_ = 0.01; // TODO: param
  }
  virtual ~ImageUpdate(){};
  int EvalRes(typename Output::Ref out, const typename Previous::CRef pre, const typename Current::CRef cur){
    for(int i=0;i<N;i++){
      if(active_[i]){
        cur.template Get<STA_BEA>()[i].Boxminus(meas_->GetBea(i),out.template Get<OUT_RES>()[i]);
      } else {
        out.template Get<OUT_RES>()[i].setZero();
      }
    }
    return 0;
  }
  int JacPre(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    return 0;
  }
  int JacCur(MatRefX J, const typename Previous::CRef pre, const typename Current::CRef cur){
    for(int i=0;i<N;i++){
      if(active_[i]){
        cur.template Get<STA_BEA>()[i].BoxminusJacInp(meas_->GetBea(i),J.block<2,2>(Output::Start(OUT_RES)+i*2,cur.Start(STA_BEA)+i*2));
      }
    }
    return 0;
  }
  virtual void AddNoise(typename Output::Ref out, MatRefX J_pre, MatRefX J_cur, const typename Previous::CRef pre, const typename Current::CRef cur){
    for(int i=0;i<N;i++){
      double w = GetWeight();
      const double norm = out.template Get<OUT_RES>()[i].norm();
      if(norm > huberTh_){
        w *= sqrt(2 * huberTh_ * (norm - 0.5 * huberTh_)/(norm*norm));
      }
      out.template Get<OUT_RES>()[i] *= w;
      J_cur.block(Output::Start(OUT_RES)+i*2,0,2,J_cur.cols()) *= w;
    }
  }
  std::array<bool,N> active_;
  double huberTh_;
};

} // namespace tsif

#endif  // TSIF_IMAGE_UPDATE_H_
