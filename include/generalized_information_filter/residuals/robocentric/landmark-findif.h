#ifndef GIF_LANDMARKFINDIF_HPP_
#define GIF_LANDMARKFINDIF_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/binary-residual.h"

namespace GIF {

/*! \brief Robocentric landmark residual.
 */
template<int NumLandmark>
class LandmarkFindif
      : public BinaryResidual<ElementPack<std::array<Vec3,NumLandmark>>,
                              ElementPack<std::array<Vec3,NumLandmark>, Vec3, Vec3>,
                              ElementPack<std::array<Vec3,NumLandmark>>,
                              ElementPack<std::array<Vec3,NumLandmark>>, EmptyMeas> {
 public:
  using mtResidual = BinaryResidual<ElementPack<std::array<Vec3,NumLandmark>>,
                                    ElementPack<std::array<Vec3,NumLandmark>, Vec3, Vec3>,
                                    ElementPack<std::array<Vec3,NumLandmark>>,
                                    ElementPack<std::array<Vec3,NumLandmark>>, EmptyMeas>;
  using mtResidual::dt_;
  using mtResidual::GetJacBlockPre;
  using mtResidual::GetJacBlockCur;
  using mtResidual::GetJacBlockNoi;
  LandmarkFindif(const std::string& name,
            const std::array<std::string,1>& innName = {"MrML"},
            const std::array<std::string,3>& preName = {"MrML", "MvM", "MwM"},
            const std::array<std::string,1>& curName = {"MrML"},
            const std::array<std::string,1>& noiName = {"MrML"})
    : mtResidual(name,innName,preName,curName,noiName,false,true,true){
    for(int i=0;i<NumLandmark;i++){
      propagation_flags_[i] = true;
    }
  }
  virtual ~LandmarkFindif() {
  }
  enum ElementsInn {FEA_INN};
  enum ElementsPre {FEA_PRE, VEL_PRE, ROR_PRE};
  enum ElementsCur {FEA_CUR};
  enum ElementsNoi {FEA_NOI};
  void Eval(std::array<Vec3,NumLandmark>& MrML_inn, const std::array<Vec3,NumLandmark>& MrML_pre,
            const Vec3& MvM_pre, const Vec3& MwM_pre,
            const std::array<Vec3,NumLandmark>& MrML_cur,
            const std::array<Vec3,NumLandmark>& MrML_noi) const {
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        MrML_inn[i] = (Mat3::Identity() - gSM(dt_ * MwM_pre)) * MrML_pre[i]
                      - dt_ * MvM_pre + MrML_noi[i] * sqrt(dt_) - MrML_cur[i];
      } else {
        MrML_inn[i] = MrML_noi[i] * sqrt(dt_);
      }
    }
  }
  void JacPre(MatX& J, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& MvM_pre, const Vec3& MwM_pre,
              const std::array<Vec3,NumLandmark>& MrML_cur,
              const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockPre<FEA_INN, FEA_PRE>(J).template block<3,3>(i*3,i*3) =
            (Mat3::Identity() - gSM(dt_ * MwM_pre));
        this->template GetJacBlockPre<FEA_INN, VEL_PRE>(J).template block<3,3>(i*3,0) =
            -dt_ * Mat3::Identity();
        this->template GetJacBlockPre<FEA_INN, ROR_PRE>(J).template block<3,3>(i*3,0) =
            dt_ * gSM(MrML_pre[i]);
      }
    }
  }
  void JacCur(MatX& J, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& MvM_pre, const Vec3& MwM_pre,
              const std::array<Vec3,NumLandmark>& MrML_cur,
              const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    for(int i=0;i<NumLandmark;i++){
      if(propagation_flags_[i]){
        this->template GetJacBlockCur<FEA_INN, FEA_CUR>(J).template block<3,3>(i*3,i*3) =
            -Mat3::Identity();
      }
    }
  }
  void JacNoi(MatX& J, const std::array<Vec3,NumLandmark>& MrML_pre,
              const Vec3& MvM_pre, const Vec3& MwM_pre,
              const std::array<Vec3,NumLandmark>& MrML_cur,
              const std::array<Vec3,NumLandmark>& MrML_noi) const {
    J.setZero();
    for(int i=0;i<NumLandmark;i++){
      this->template GetJacBlockNoi<FEA_INN, FEA_NOI>(J).template block<3,3>(i*3,i*3) =
          sqrt(dt_) * Mat3::Identity();
    }
  }
  void SetPropagationFlag(int i, bool flag) {
    propagation_flags_[i] = flag;
  }

 protected:
  std::array<bool,NumLandmark> propagation_flags_;
};

}
#endif /* GIF_LANDMARKFINDIF_HPP_ */
