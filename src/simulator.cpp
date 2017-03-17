#include "tsif/utils/simulator.h"

Simulator::Simulator(){
  // Simulated parameters
  sim_bw_ = tsif::Vec3(0.2,-0.1,0.1);
  sim_bf_ = tsif::Vec3(0.3,0.1,-0.1);
  Ig_ = tsif::Vec3(0,0,-9.81);
  sim_IrIJ_ = tsif::Vec3(0,0,0);
  sim_qIJ_ = tsif::Quat(1,0,0,0);
  sim_MrMC_ = tsif::Vec3(-0.0111674199187,-0.0574640920022,0.0207586947896);
  sim_qMC_ = tsif::Quat(-0.712115587266,0.00666398307551,-0.0079168224269,-0.701985972528);
  sim_qIM_des_ = tsif::Quat(1,0,0,0);

  // Simulation parameters
  sim_seed_ = 0;
  sim_dt_ = 0.01;
  sim_aa_ = 1;
  sim_ar_ = 0.1;
  sim_av_ = 1;
  sim_an_ = 5;
  sim_wq_ = 0.1;
  sim_ww_ = 0.1;
  sim_wn_ = 0.5;

  sim_noise_gyr_ = 2e-3;
  sim_outlier_noise_gyr_ = 1;
  sim_outlier_rate_gyr_ = 1000;
  sim_noise_acc_ = 1e-2;
  sim_outlier_noise_acc_ = 10;
  sim_outlier_rate_acc_ = 1000;
  sim_noise_pos_ = 1e-2;
  sim_outlier_noise_pos_ = 10;
  sim_outlier_rate_pos_ = 1000;
  sim_noise_rot_ = 1e-2;
  sim_landmark_spread_ = 10;
  allowOutlier_ = false;

  for(int i=0;i<kL;i++){
    JrJL_[i] = sim_landmark_spread_*tsif::NormalRandomNumberGenerator::Instance().GetVec<3>();
  }

  init();
}

Simulator::~Simulator(){
}

void Simulator::step(){
  tsif::NormalRandomNumberGenerator& gen = tsif::NormalRandomNumberGenerator::Instance();
  t_ = t_ + sim_dt_;
  sim_IaM_ = (1-sim_aa_*sim_dt_) * sim_IaM_
             - sim_ar_*sim_dt_ * sim_IrIM_
             - sim_av_*sim_dt_ * sim_IvM_
             + sim_an_*std::sqrt(sim_dt_)*gen.GetVec<3>();
  sim_IvM_ = sim_IvM_ + sim_dt_ * sim_IaM_;
  sim_IrIM_ = sim_IrIM_ + sim_dt_ * sim_IvM_;

  sim_IwIM_ = (1-sim_ww_*sim_dt_) * sim_IwIM_
             + sim_wq_*sim_dt_ * tsif::Boxminus(sim_qIM_des_,sim_qIM_)
             + sim_wn_*std::sqrt(sim_dt_)*gen.GetVec<3>();
  tsif::Quat sim_qIM_next = tsif::Boxplus(sim_qIM_,sim_IwIM_*sim_dt_);
  sim_qIM_ = sim_qIM_next;

  genMeas();
}
void Simulator::init(){
  tsif::NormalRandomNumberGenerator::Instance().SetSeed(sim_seed_);

  // Simulated trajectory
  t_ = 0;
  sim_IrIM_ = tsif::Vec3(0,0,0);
  sim_qIM_ = tsif::Quat(1,0,0,0);
  sim_IvM_ = tsif::Vec3(0,0,0);
  sim_IwIM_ = tsif::Vec3(0,0,0);
  sim_IaM_ = tsif::Vec3(0,0,0);

  // Generate measurements
  genMeas();
}
void Simulator::genMeas(){
  tsif::NormalRandomNumberGenerator& gen = tsif::NormalRandomNumberGenerator::Instance();

  // IMU measurement
  std::bernoulli_distribution distribution_gyr(1/sim_outlier_rate_gyr_);
  double noise_gyr = allowOutlier_ && distribution_gyr(gen.GetGenerator()) ? sim_outlier_noise_gyr_ : sim_noise_gyr_;
  meas_MwM_ = sim_qIM_.toRotationMatrix().transpose()*(sim_IwIM_) + sim_bw_ + noise_gyr/std::sqrt(sim_dt_)*gen.GetVec<3>();
  std::bernoulli_distribution distribution_acc(1/sim_outlier_rate_acc_);
  double noise_acc = allowOutlier_ && distribution_acc(gen.GetGenerator()) ? sim_outlier_noise_acc_ : sim_noise_acc_;
  meas_MfM_ = sim_qIM_.toRotationMatrix().transpose()*(tsif::Vec3(sim_IaM_ - Ig_)) + sim_bf_ + noise_acc/std::sqrt(sim_dt_)*gen.GetVec<3>();

  // IMU measurement
  noise_gyr = allowOutlier_ && distribution_gyr(gen.GetGenerator()) ? sim_outlier_noise_gyr_ : sim_noise_gyr_;
  meas_MwM2_ = sim_qIM_.toRotationMatrix().transpose()*(sim_IwIM_) + sim_bw_ + noise_gyr/std::sqrt(sim_dt_)*gen.GetVec<3>();
  noise_acc = allowOutlier_ && distribution_acc(gen.GetGenerator()) ? sim_outlier_noise_acc_ : sim_noise_acc_;
  meas_MfM2_ = sim_qIM_.toRotationMatrix().transpose()*(tsif::Vec3(sim_IaM_ - Ig_)) + sim_bf_ + noise_acc/std::sqrt(sim_dt_)*gen.GetVec<3>();

  // Pose measurement
  std::bernoulli_distribution distribution_pos(1/sim_outlier_rate_pos_);
  double noise_pos = allowOutlier_ && distribution_pos(gen.GetGenerator()) ? sim_outlier_noise_pos_ : sim_noise_pos_;
  const tsif::Vec3 JrJC = sim_qIJ_.toRotationMatrix().transpose()*(tsif::Vec3(sim_IrIM_ - sim_IrIJ_ + sim_qIM_.toRotationMatrix()*(sim_MrMC_)));
  meas_JrJC_ = JrJC+noise_pos*gen.GetVec<3>();
  meas_qJC_ = tsif::Boxplus((sim_qIJ_.inverse()*sim_qIM_*sim_qMC_),sim_noise_rot_*gen.GetVec<3>());

  // Image measurement
  for(int i=0;i<kL;i++){
    const tsif::Vec3 CrCL = (sim_qMC_.inverse()*sim_qIM_.inverse()*sim_qIJ_).toRotationMatrix()*(JrJL_[i]-JrJC);
    meas_CcCL_isVisible_[i] = cam_.BearingToPixel(CrCL,meas_CcCL_[i]);
    meas_CdCL_[i] = CrCL.norm();
  }
}
