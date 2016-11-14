#ifndef GIF_LEGKINEMATICUPDATE_HPP_
#define GIF_LEGKINEMATICUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"

namespace GIF {

/*! \brief Leg Kinematic Measurement
 *         ElementVector that can be used to hold leg kinematic measurements. It is composed of an
 *         array of vectors containing the single encoder measurements.
 */
template<int NumLeg, int NumDof>
class KinematicMeasurement : public ElementVector {
 public:
  typedef std::array<Vec<NumDof>,NumLeg> Kin;
  KinematicMeasurement(const Kin& kin = ElementTraits<Kin>::Identity())
      : ElementVector(std::shared_ptr<ElementVectorDefinition>(new ElementPack<Kin>({"kin"}))),
        kin_(ElementVector::GetValue<Kin>("kin")) {
    kin_ = kin;
    for(int i=0;i<NumLeg;i++){
      contact_flag_[i] = false;
    }
  }
  Kin& kin_;
  std::array<bool,NumLeg> contact_flag_;
};

/*! \brief Leg Kinematic Update
 *         Transforms the kinematic measurements into relative measurement between main body and
 *         foothold. Here the foothold is directly stored in robocentric coordinates and therefore
 *         the residual condenses down to prediction - forward kinematics. The forward kinematics
 *         are expressed with respect to B (body frame) and the state is expressed with respect to
 *         M (Imu).
 *
 *         Coordinate frames:
 *           B: Body
 *           M: IMU
 *           F: Foothold
 */
template<typename KinModel>
class LegKinematicUpdate : public UnaryUpdate<ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                              ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                              ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                              KinematicMeasurement<KinModel::kNumLeg,
                                                                   KinModel::kNumDof>> {
 public:
  using mtUnaryUpdate = UnaryUpdate<ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                    ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                    ElementPack<std::array<Vec3,KinModel::kNumLeg>>,
                                    KinematicMeasurement<KinModel::kNumLeg,
                                                         KinModel::kNumDof>>;
  using mtUnaryUpdate::meas_;
  typedef std::array<Vec3,KinModel::kNumLeg> LegArray;
  enum Elements {KIN};

  LegKinematicUpdate(const std::string& name,
                     const std::string& errorName = "BrBF",
                     const std::string& stateName = "MrMF",
                     const std::string& noiseName = "BrBF")
      : mtUnaryUpdate(name,{errorName},{stateName},{noiseName}) {
  }

  virtual ~LegKinematicUpdate() {
  }

  void Eval(LegArray& MrMF_inn, const LegArray& MrMF_cur, const LegArray& MrMF_noi) const {
    for(int i=0;i<KinModel::kNumLeg;i++){
      if(GetContactFlagFromMeas(i)){
        const Vec3 BrBF = BrBM_ + qMB_.inverseRotate(MrMF_cur[i]);
        MrMF_inn[i] = BrBF
            - model_->forwardKinematicsBaseToFootInBaseFrame(meas_->kin_[i],i) + MrMF_noi[i];
      } else {
        MrMF_inn[i] = MrMF_noi[i];
      }
    }
  }
  void JacCur(MatX& J, const LegArray& MrMF_cur, const LegArray& MrMF_noi) const {
    J.setZero();
    for(int i=0;i<KinModel::kNumLeg;i++){
      if(GetContactFlagFromMeas(i)){
        this->template GetJacBlockCur<KIN, KIN>(J).template block<3,3>(3*i,3*i) =
            RotMat(qMB_.inverted()).matrix();
      }
    }
  }
  void JacNoi(MatX& J, const LegArray& MrMF_cur, const LegArray& MrMF_noi) const {
    J.setZero();
    for(int i=0;i<KinModel::kNumLeg;i++){
      this->template GetJacBlockNoi<KIN, KIN>(J).template block<3,3>(3*i,3*i) = Mat3::Identity();
    }
  }

  void SetModelPtr(const std::shared_ptr<KinModel>& model) {
    model_ = model;
  }

  void SetExtrinsics(Vec3 BrBM, Quat qMB){
    BrBM_  = BrBM;
    qMB_ = qMB;
  }

  bool GetContactFlagFromMeas(int i) const{
    return meas_->contact_flag_[i];
  }

  Vec3 GetMrMFFromMeas(int i) const{
    return qMB_.rotate(Vec3(model_->forwardKinematicsBaseToFootInBaseFrame(meas_->kin_[i],i)
                            - BrBM_));
  }

 protected:
  std::shared_ptr<KinModel> model_;
  Vec3 BrBM_;
  Quat qMB_;
};

}

#endif /* GIF_LEGKINEMATICUPDATE_HPP_ */
