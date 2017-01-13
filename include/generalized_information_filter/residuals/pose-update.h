#ifndef GIF_POSEUPDATE_HPP_
#define GIF_POSEUPDATE_HPP_

#include "generalized_information_filter/common.h"
#include "generalized_information_filter/unary-update.h"
#include "generalized_information_filter/measurements/pose-meas.h"

namespace GIF {

template<int i, typename ... Ts>
struct MultiPosePack{
  typedef typename MultiPosePack<i-1,Vec3,Quat,Ts...>::type type;
};

template<typename ... Ts>
struct MultiPosePack<0,Ts...>{
  typedef ElementPack<Ts...> type;
};

/*! \brief Pose Update
 *         Builds a residual between current estimated pose and an external measured pose. The state
 *         is parametrized by the position (IrIB) and attitude (qIB). Optionally, the offset
 *         between the inertial frames (IrIJ and qIJ) and the extrinsic calibration of the pose
 *         measurement (BrBC and qBC) are also co-estimated.
 *
 *         Coordinate frames:
 *           B: Body
 *           C: Pose measurement sensor
 *           J: Pose measurement reference inertial frame
 *           I: Estimation reference inertial frame
 */
template<bool doInertialAlignment, bool doBodyAlignment>
class PoseUpdate : public UnaryUpdate<ElementPack<Vec3, Quat>,
    typename MultiPosePack<1+doInertialAlignment+doBodyAlignment>::type,
    ElementPack<Vec3, Vec3>, PoseMeas> {
 public:
  static constexpr int numPose = 1+doInertialAlignment+doBodyAlignment;
  static constexpr int POS_INN = 0;
  static constexpr int ATT_INN = 1;
  static constexpr int POS_CUR = 0;
  static constexpr int ATT_CUR = 1;
  static constexpr int IJP_CUR = POS_CUR+2*doInertialAlignment;
  static constexpr int IJA_CUR = ATT_CUR+2*doInertialAlignment;
  static constexpr int BCP_CUR = IJP_CUR+2*doBodyAlignment;
  static constexpr int BCA_CUR = IJA_CUR+2*doBodyAlignment;
  static constexpr int POS_NOI = 0;
  static constexpr int ATT_NOI = 1;
  using mtUnaryUpdate = UnaryUpdate<ElementPack<Vec3, Quat>,
      typename MultiPosePack<1+doInertialAlignment+doBodyAlignment>::type,
      ElementPack<Vec3, Vec3>, PoseMeas>;
  using mtUnaryUpdate::meas_;

  PoseUpdate(const std::string& name,
             const std::array<std::string,2>& errorName,
             const std::array<std::string,2*numPose>& stateName,
             const std::array<std::string,2>& noiseName)
       : mtUnaryUpdate(name, errorName, stateName, noiseName){
    BrBC_.setZero();
    qBC_.setIdentity();
    IrIJ_.setZero();
    qIJ_.setIdentity();
    useAttitude_ = true;
    usePosition_ = true;
    huberTh_ = -1.0;
  }

  virtual ~PoseUpdate() {
  }

  // Full version
  void Eval(Vec3& JrJC_inn, Quat& qJC_inn,
            const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
            const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    if(usePosition_){
      JrJC_inn = qIJ_cur.inverseRotate(Vec3(IrIB_cur - IrIJ_cur + qIB_cur.rotate(BrBC_cur)))
          - meas_->JrJC_ + JrJC_noi;
    } else {
      JrJC_inn = JrJC_noi;
    }
    Quat dQ = dQ.exponentialMap(qJC_noi);
    if(useAttitude_){
      qJC_inn = dQ * qIJ_cur.inverted() * qIB_cur * qBC_cur * meas_->qJC_.inverted();
    } else {
      qJC_inn = dQ;
    }
  }

  // State only version
  void Eval(Vec3& JrJC_inn, Quat& qJC_inn, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    Eval(JrJC_inn, qJC_inn, IrIB_cur, qIB_cur, IrIJ_, qIJ_, BrBC_, qBC_, JrJC_noi, qJC_noi);
  }

  // Single calibration version
  void Eval(Vec3& JrJC_inn, Quat& qJC_inn, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& r, const Quat& q, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    if(doInertialAlignment){
      Eval(JrJC_inn, qJC_inn, IrIB_cur, qIB_cur, r, q, BrBC_, qBC_, JrJC_noi, qJC_noi);
    } else {
      Eval(JrJC_inn, qJC_inn, IrIB_cur, qIB_cur, IrIJ_, qIJ_, r, q, JrJC_noi, qJC_noi);
    }
  }

  // Full version
  void JacCur(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    if(usePosition_){
      this->template GetJacBlockCur<POS_INN, POS_CUR>(J) = RotMat(qIJ_cur).matrix().transpose();
      this->template GetJacBlockCur<POS_INN, ATT_CUR>(J) = -gSM(RotMat(qIJ_cur.inverted()
        * qIB_cur).rotate(BrBC_cur)) * RotMat(qIJ_cur).matrix().transpose();
      if(doInertialAlignment){
        this->template GetJacBlockCur<POS_INN, IJP_CUR>(J) = -RotMat(qIJ_cur).matrix().transpose();
        this->template GetJacBlockCur<POS_INN, IJA_CUR>(J) = RotMat(qIJ_cur).matrix().transpose()
          * gSM(Vec3(IrIB_cur - IrIJ_cur + qIB_cur.rotate(BrBC_cur)));
      }
      if(doBodyAlignment){
        this->template GetJacBlockCur<POS_INN, BCP_CUR>(J) = RotMat(qIJ_cur.inverted()*qIB_cur).matrix();
      }
    }
    if(useAttitude_){
      this->template GetJacBlockCur<ATT_INN, ATT_CUR>(J) = RotMat(qIJ_cur.inverted()).matrix();
      if(doInertialAlignment){
        this->template GetJacBlockCur<ATT_INN, IJA_CUR>(J) = -RotMat(qIJ_cur.inverted()).matrix();
      }
      if(doBodyAlignment){
        this->template GetJacBlockCur<ATT_INN, BCA_CUR>(J) = RotMat(qIJ_cur.inverted() * qIB_cur).matrix();
      }
    }
  }

  // State only version
  void JacCur(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    JacCur(J, IrIB_cur, qIB_cur, IrIJ_, qIJ_, BrBC_, qBC_, JrJC_noi, qJC_noi);
  }

  // Single calibration version
  void JacCur(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& r, const Quat& q, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    if(doInertialAlignment){
      JacCur(J, IrIB_cur, qIB_cur, r, q, BrBC_, qBC_, JrJC_noi, qJC_noi);
    } else {
      JacCur(J, IrIB_cur, qIB_cur, IrIJ_, qIJ_, r, q, JrJC_noi, qJC_noi);
    }
  }

  // Full version
  void JacNoi(MatX& J,
              const Vec3& IrIB_cur, const Quat& qIB_cur, const Vec3& IrIJ_cur, const Quat& qIJ_cur,
              const Vec3& BrBC_cur, const Quat& qBC_cur, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    J.setZero();
    this->template GetJacBlockNoi<POS_INN, POS_NOI>(J) = Mat3::Identity();
    this->template GetJacBlockNoi<ATT_INN, ATT_NOI>(J) = Mat3::Identity();
  }

  // State only version
  void JacNoi(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    JacNoi(J, IrIB_cur, qIB_cur, IrIJ_, qIJ_, BrBC_, qBC_, JrJC_noi, qJC_noi);
  }

  // Single calibration version
  void JacNoi(MatX& J, const Vec3& IrIB_cur, const Quat& qIB_cur,
            const Vec3& r, const Quat& q, const Vec3& JrJC_noi, const Vec3& qJC_noi) const {
    if(doInertialAlignment){
      JacNoi(J, IrIB_cur, qIB_cur, r, q, BrBC_, qBC_, JrJC_noi, qJC_noi);
    } else {
      JacNoi(J, IrIB_cur, qIB_cur, IrIJ_, qIJ_, r, q, JrJC_noi, qJC_noi);
    }
  }

  void SetInertialAlignment(const Vec3& IrIJ, const Quat& qIJ){
    IrIJ_  = IrIJ;
    qIJ_ = qIJ;
  }
  void SetBodyAlignment(const Vec3& BrBC, const Quat& qBC){
    BrBC_  = BrBC;
    qBC_ = qBC;
  }
  void SetPositionFlag(bool mb){
    usePosition_ = mb;
  }
  void SetAttitudeFlag(bool mb){
    useAttitude_ = mb;
  }
  void SetHuberTh(double th){
    huberTh_ = th;
  }
  double GetNoiseWeighting(const ElementVector& inn, int i){
    if(huberTh_ >= 0){
      double norm = inn.GetValue<Vec3>(POS_INN).norm(); // TODO: position only so far
      if(norm > huberTh_){
        LOG(WARNING) << "Outlier on position update: " << norm << std::endl;
        return sqrt(huberTh_ * (norm - 0.5 * huberTh_)/(norm*norm));
      } else {
        return 1.0;
      }
    } else {
      return 1.0;
    }
  }

 protected:
  Vec3 BrBC_;
  Quat qBC_;
  Vec3 IrIJ_;
  Quat qIJ_;
  bool useAttitude_;
  bool usePosition_;
  double huberTh_;
};

}

#endif /* GIF_POSEUPDATE_HPP_ */
