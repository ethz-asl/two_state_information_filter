#include "tsif/utils/simulator.h"
#include "tsif/filters/vio.h"

int main(int argc, char** argv){
  std::string optionFile = "/home/michael/workspace/generalized_information_filter/cfg/vio.cfg";
  std::string taskFile = "/home/michael/workspace/generalized_information_filter/cfg/vio_task.cfg";
  bool useSim = tsif::OptionLoader::Instance().Get<int>(optionFile,"use_sim");
  const int L = 25;
  tsif::VioFilter<L> filter(optionFile);
  filter.GetState().SetIdentity();
  Simulator sim;
  sim.allowOutlier_ = false;
  sim.init();

  std::string folder = tsif::OptionLoader::Instance().Get<std::string>(taskFile,"data_folder");

  std::ifstream cam0data(folder + "/mav0/cam0/data.csv");
  std::string cam0string;
  getline(cam0data, cam0string); // Ignore first
  double cam0time;

  std::ifstream imu0data(folder + "/mav0/imu0/data.csv");
  std::string imu0string;
  getline(imu0data, imu0string); // Ignore first
  double imu0time = std::numeric_limits<double>::min();

  for(int i=0;i<10000 && cam0data.good();i++){
    if(useSim){
      sim.step();
      std::shared_ptr<tsif::MeasImg<L>> meas(new tsif::MeasImg<L>());
      for(int i=0;i<Simulator::kL;i++){
        meas->isSim_ = true;
        if(sim.meas_CcCL_isVisible_[i] && sim.meas_CcCL_[i](0) >= 0 && sim.meas_CcCL_[i](0) < 752 && sim.meas_CcCL_[i](1) >= 0 && sim.meas_CcCL_[i](1) < 480){
          meas->keyPoints_.push_back(cv::KeyPoint(sim.meas_CcCL_[i](0),sim.meas_CcCL_[i](1),0,sim.meas_CdCL_[i],0,0,i));
        }
      }
      filter.AddMeas<2>(tsif::TimePoint(tsif::fromSec(sim.t_+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetAcc"))),std::make_shared<tsif::MeasAcc>(sim.meas_MfM_));
      filter.AddMeas<3>(tsif::TimePoint(tsif::fromSec(sim.t_+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetGyr"))),std::make_shared<tsif::MeasGyr>(sim.meas_MwM_));
      filter.AddMeas<6>(tsif::TimePoint(tsif::fromSec(sim.t_+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetCam"))),meas);
    } else {
      getline (cam0data, cam0string, ',');
      if(!cam0data.good()){
        break;
      }
      cam0time = stod(cam0string)*1e-9;
      getline (cam0data, cam0string);
      if(cam0string.back() == '\r'){
        cam0string.pop_back();
      }
      std::string imgFilename = folder + "/mav0/cam0/data/" + cam0string;
      std::shared_ptr<tsif::MeasImg<L>> meas(new tsif::MeasImg<L>(imgFilename.c_str()));
      for(int i=0;i<Simulator::kL;i++){
        if(sim.meas_CcCL_isVisible_[i] && sim.meas_CcCL_[i](0) >= 0 && sim.meas_CcCL_[i](0) < 752 && sim.meas_CcCL_[i](1) >= 0 && sim.meas_CcCL_[i](1) < 480){
          meas->keyPoints_.push_back(cv::KeyPoint(sim.meas_CcCL_[i](0),sim.meas_CcCL_[i](1),0,sim.meas_CdCL_[i],0,0,i));
        }
      }
      filter.AddMeas<6>(tsif::TimePoint(tsif::fromSec(cam0time+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetCam"))),meas);

      while(imu0time < cam0time){
        getline (imu0data, imu0string, ',');
        if(!imu0data.good()){
          break;
        }
        imu0time = stod(imu0string)*1e-9;
        getline (imu0data, imu0string, ',');
        const double wx = stod(imu0string);
        getline (imu0data, imu0string, ',');
        const double wy = stod(imu0string);
        getline (imu0data, imu0string, ',');
        const double wz = stod(imu0string);
        getline (imu0data, imu0string, ',');
        const double ax = stod(imu0string);
        getline (imu0data, imu0string, ',');
        const double ay = stod(imu0string);
        getline (imu0data, imu0string);
        const double az = stod(imu0string);
        std::shared_ptr<tsif::MeasGyr> measGyr(new tsif::MeasGyr(tsif::Vec3(wx,wy,wz)));
        std::shared_ptr<tsif::MeasAcc> measAcc(new tsif::MeasAcc(tsif::Vec3(ax,ay,az)));
        filter.AddMeas<2>(tsif::TimePoint(tsif::fromSec(imu0time+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetAcc"))),measAcc);
        filter.AddMeas<3>(tsif::TimePoint(tsif::fromSec(imu0time+tsif::OptionLoader::Instance().Get<double>(optionFile,"timeOffsetGyr"))),measGyr); // TODO timing
      }
    }

    filter.Update();
//    std::cout << "====================\n" << filter.GetState().Print();

    //    const tsif::Vec3 IrIM = sim.sim_IrIM_;
    //    const tsif::Quat qIM = sim.sim_qIM_;
    //    const tsif::Vec3 vel = qIM.toRotationMatrix().transpose()*sim.sim_IvM_;
    //    const tsif::Vec3 ror = qIM.toRotationMatrix().transpose()*sim.sim_IwIM_;
    //    const tsif::Vec3 JrJC = sim.sim_qIJ_.toRotationMatrix().transpose()*(tsif::Vec3(IrIM - sim.sim_IrIJ_ + qIM.toRotationMatrix()*(sim.sim_MrMC_)));
    //    const tsif::Vec3 CrCL = (sim.sim_qMC_.inverse()*qIM.inverse()*sim.sim_qIJ_).toRotationMatrix()*(sim.JrJL_[i]-JrJC);
    //    const double invDis = 1/CrCL.norm();
    //    const tsif::UnitVector CnCL(CrCL);
    //
    //    const double d = 1e-8;
    //    const tsif::Vec3 IrIM_pert = IrIM + d*qIM.toRotationMatrix()*vel;
    //    const tsif::Quat qIM_pert = qIM * tsif::Exp(ror*d);
    //    const tsif::Vec3 JrJC_pert = sim.sim_qIJ_.toRotationMatrix().transpose()*(tsif::Vec3(IrIM_pert - sim.sim_IrIJ_ + qIM_pert.toRotationMatrix()*(sim.sim_MrMC_)));
    //    const tsif::Vec3 CrCL_pert = (sim.sim_qMC_.inverse()*qIM_pert.inverse()*sim.sim_qIJ_).toRotationMatrix()*(sim.JrJL_[i]-JrJC_pert);
    //    const double invDis_pert = 1/CrCL_pert.norm();
    //    const tsif::UnitVector CnCL_pert(CrCL_pert);
    //
    //    tsif::Vec<2> dif;
    //    CnCL_pert.Boxminus(CnCL,dif);
    //    std::cout << dif.transpose()/d << std::endl;
    //    std::cout << (CnCL.GetN().transpose()*(ror + invDis * CnCL.GetVec().cross(vel))).transpose() << std::endl;
    //    std::cout << (invDis_pert-invDis)/d << std::endl;
    //    std::cout << CnCL.GetVec().dot(vel)*invDis*invDis << std::endl;
  }
  std::cout << filter.PrintConnectivity();
  filter.JacTestAll(1e-6,1e-8);
  std::cout << "=================== State ===================" << std::endl;
  std::cout << filter.GetState().Print();
  std::cout << "=================== GT ===================" << std::endl;
  std::cout << sim.meas_JrJC_.transpose() << std::endl;
  std::cout << sim.meas_qJC_.coeffs().transpose() << std::endl;
  std::cout << (sim.sim_qIM_.toRotationMatrix().transpose()*sim.sim_IvM_).transpose() << std::endl;
  std::cout << (sim.sim_qIM_.toRotationMatrix().transpose()*sim.sim_IwIM_).transpose() << std::endl;

  return 0;
}
