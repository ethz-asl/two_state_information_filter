#include "tsif/utils/simulator.h"
#include "tsif/filters/vio.h"

int main(int argc, char** argv){
  tsif::VioFilter filter;
  filter.GetState().SetIdentity();
  Simulator sim;
  sim.allowOutlier_ = false;
  sim.init();

  for(int i=0;i<1;i++){
    sim.step();
    cv::Mat img = cv::Mat::zeros(480, 752, CV_8U);

    std::shared_ptr<tsif::MeasImg<10>> meas(new tsif::MeasImg<10>());
    for(int i=0;i<Simulator::kL;i++){
      if(sim.meas_CcCL_isVisible_[i]){
        meas->keyPoints_.push_back(cv::KeyPoint(sim.meas_CcCL_[i](0),sim.meas_CcCL_[i](1),0,0,0,0,i));
      }
    }
    cv::drawKeypoints(img,meas->keyPoints_,img);

    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Display Image", img);
    cv::waitKey(1);

    filter.AddMeas<2>(tsif::TimePoint(tsif::fromSec(sim.t_)),std::make_shared<tsif::MeasAcc>(sim.meas_MfM_));
    filter.AddMeas<3>(tsif::TimePoint(tsif::fromSec(sim.t_-sim.sim_dt_)),std::make_shared<tsif::MeasGyr>(sim.meas_MwM_));
    filter.AddMeas<6>(tsif::TimePoint(tsif::fromSec(sim.t_)),meas);

    filter.Update();
  }
  std::cout << "=================== State ===================" << std::endl;
  std::cout << filter.GetState().Print();
  std::cout << filter.PrintConnectivity();
  filter.JacTestAll(1e-6,1e-8);

  return 0;
}
