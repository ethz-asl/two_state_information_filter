#include "tsif/filters/imu_gps.h"
#include "tsif/utils/simulator.h"

int main(int argc, char** argv){
  tsif::ImuGpsFilter filter;
  Simulator sim;
  sim.allowOutlier_ = false;
  sim.init();

  for(int i=0;i<10000;i++){
    sim.step();
    filter.AddMeas<2>(tsif::TimePoint(tsif::fromSec(sim.t_)),std::make_shared<tsif::MeasAcc>(sim.meas_MfM_));
    filter.AddMeas<3>(tsif::TimePoint(tsif::fromSec(sim.t_-sim.sim_dt_)),std::make_shared<tsif::MeasGyr>(sim.meas_MwM_));
    filter.AddMeas<6>(tsif::TimePoint(tsif::fromSec(sim.t_)),std::make_shared<tsif::MeasPos>(sim.meas_JrJC_));
    filter.AddMeas<8>(tsif::TimePoint(tsif::fromSec(sim.t_)),std::make_shared<tsif::MeasAtt>(sim.meas_qJC_));
    filter.Update();
  }
  std::cout << "=================== State ===================" << std::endl;
  std::cout << filter.GetState().Print();
  std::cout << "=================== GT ===================" << std::endl;
  std::cout << sim.meas_JrJC_.transpose() << std::endl;
  std::cout << sim.meas_qJC_.coeffs().transpose() << std::endl;
  std::cout << (sim.sim_qIM_.toRotationMatrix().transpose()*sim.sim_IvM_).transpose() << std::endl;
  std::cout << (sim.sim_qIM_.toRotationMatrix().transpose()*sim.sim_IwIM_).transpose() << std::endl;
  std::cout << filter.PrintConnectivity();
  filter.JacTestAll(1e-6,1e-8);

  return 0;
}
