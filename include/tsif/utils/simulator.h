#ifndef TSIF_SIMULATOR_H_
#define TSIF_SIMULATOR_H_

#include "tsif/utils/common.h"
#include "tsif/utils/camera.h"

class Simulator{
public:
  // Simulated parameters
  tsif::Vec3 sim_bw_;
  tsif::Vec3 sim_bf_;
  tsif::Vec3 Ig_;
  tsif::Vec3 sim_IrIJ_;
  tsif::Quat sim_qIJ_;
  tsif::Vec3 sim_MrMC_;
  tsif::Quat sim_qMC_;
  tsif::Quat sim_qIM_des_;

  // Simulated trajectory
  double t_;
  tsif::Vec3 sim_IrIM_;
  tsif::Quat sim_qIM_;
  tsif::Vec3 sim_IvM_;
  tsif::Vec3 sim_IwIM_;
  tsif::Vec3 sim_IaM_;

  // IMU measurement
  tsif::Vec3 meas_MwM_;
  tsif::Vec3 meas_MfM_;
  tsif::Vec3 meas_MwM2_;
  tsif::Vec3 meas_MfM2_;

  // Pose measurement
  tsif::Vec3 meas_JrJC_;
  tsif::Quat meas_qJC_;

  // Landmarks
  tsif::Camera cam_;
  static const int kL = 1000;
  std::array<tsif::Vec3,kL> JrJL_;
  std::array<tsif::Vec<2>,kL> meas_CcCL_;
  std::array<double,kL> meas_CdCL_;
  std::array<bool,kL> meas_CcCL_isVisible_;

  Simulator();
  ~Simulator();
  void step();
  void init();
  void genMeas();

  // Simulation parameters
  int sim_seed_;
  double sim_dt_;
  double sim_aa_;
  double sim_ar_;
  double sim_av_;
  double sim_an_;
  double sim_ww_;
  double sim_wq_;
  double sim_wn_;
  double sim_noise_gyr_;
  double sim_outlier_noise_gyr_;
  double sim_outlier_rate_gyr_;
  double sim_noise_acc_;
  double sim_outlier_noise_acc_;
  double sim_outlier_rate_acc_;
  double sim_noise_pos_;
  double sim_outlier_noise_pos_;
  double sim_outlier_rate_pos_;
  double sim_noise_rot_;
  double sim_landmark_spread_;
  bool allowOutlier_;
};

#endif /* TSIF_SIMULATOR_H_ */
