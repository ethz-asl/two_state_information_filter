#include "tsif/utils/simulator.h"
#include "tsif/filters/vio.h"

int main(int argc, char** argv){
  const int L = 10;
  tsif::VioFilter<L> filter;
  filter.GetState().SetIdentity();
  Simulator sim;
  sim.allowOutlier_ = false;
  sim.init();

  for(int i=0;i<1000;i++){
    sim.step();
    cv::Mat img = cv::Mat::zeros(480, 752, CV_8U);


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

    std::shared_ptr<tsif::MeasImg<L>> meas(new tsif::MeasImg<L>());
    for(int i=0;i<Simulator::kL;i++){
      if(sim.meas_CcCL_isVisible_[i] && sim.meas_CcCL_[i](0) >= 0 && sim.meas_CcCL_[i](0) < 752 && sim.meas_CcCL_[i](1) >= 0 && sim.meas_CcCL_[i](1) < 480){
        meas->keyPoints_.push_back(cv::KeyPoint(sim.meas_CcCL_[i](0),sim.meas_CcCL_[i](1),0,sim.meas_CdCL_[i],0,0,i));
      }
    }
//    cv::drawKeypoints(img,meas->keyPoints_,img);

//    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE);
//    cv::imshow("Display Image", img);
//    cv::waitKey(0);

    filter.AddMeas<2>(tsif::TimePoint(tsif::fromSec(sim.t_)),std::make_shared<tsif::MeasAcc>(sim.meas_MfM_));
    filter.AddMeas<3>(tsif::TimePoint(tsif::fromSec(sim.t_-sim.sim_dt_)),std::make_shared<tsif::MeasGyr>(sim.meas_MwM_));
    filter.AddMeas<6>(tsif::TimePoint(tsif::fromSec(sim.t_)),meas);

    filter.Update();
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
