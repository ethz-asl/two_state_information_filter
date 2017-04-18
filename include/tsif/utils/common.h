#ifndef TSIF_COMMON_HPP_
#define TSIF_COMMON_HPP_

#include <array>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <fstream>

#if TSIF_VERBOSE > 0
#define TSIF_LOG(msg) std::cout << msg << std::endl
#define TSIF_LOGIF(con,msg) if(con) TSIF_LOG(msg)
#define TSIF_LOGW(msg) std::cout << "\033[33m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGWIF(con,msg) if(con) TSIF_LOGW(msg)
#define TSIF_LOGE(msg) std::cout << "\033[31m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGEIF(con,msg) if(con) TSIF_LOGE(msg)
#else
#define TSIF_LOG(msg)
#define TSIF_LOGIF(con,msg)
#define TSIF_LOGW(msg)
#define TSIF_LOGWIF(con,msg)
#define TSIF_LOGE(msg)
#define TSIF_LOGEIF(con,msg)
#endif

namespace tsif{

template<int N = -1>
using Vec = Eigen::Matrix<double,N,1>;
using Vec3 = Vec<3>;
using VecX = Vec<>;
template<int N = -1>
using VecRef = Eigen::Ref<Vec<N>>;
using VecRef3 = VecRef<3>;
using VecRefX = VecRef<>;
template<int N = -1>
using VecCRef = Eigen::Ref<const Vec<N>>;
using VecCRef3 = VecCRef<3>;
using VecCRefX = VecCRef<>;

template<int N = -1, int M = N>
using Mat = Eigen::Matrix<double,N,M>;
using Mat3 = Mat<3>;
using MatX = Mat<>;
template<int N = -1, int M = N>
using MatRef = typename std::conditional<(N==1 & M>1),
                                         Eigen::Ref<Mat<N,M>,0,Eigen::InnerStride<>>,
                                         Eigen::Ref<Mat<N,M>>>::type;
using MatRef3 = MatRef<3>;
using MatRefX = MatRef<>;
template<int N = -1, int M = N>
using MatCRef = Eigen::Ref<const Mat<N,M>>;
using MatCRef3 = MatCRef<3>;
using MatCRefX = MatCRef<>;

typedef Eigen::Quaterniond Quat;
static constexpr double Sinc(const double x){
  return fabs(x) < 1e-8 ? 1 : sin(x)/x;
}
static Quat Exp(const Vec3& v){
  const double ha = 0.5*v.norm();
  const double re = cos(ha);
  const Vec3 im = 0.5*Sinc(ha)*v;
  return Quat(re,im(0),im(1),im(2));
}
static Vec3 Log(const Quat& q){
  const double re = q.w();
  const Vec3 im(q.x(),q.y(),q.z());
  const double sha = im.norm();
  return sha < 1e-8 ? ((std::signbit(re)!=0)*-2+1)*2*im : 2*atan2(sha,re)/sha*im;
}
static Quat Boxplus(const Quat& q, const Vec3& v){
  return Exp(v)*q;
}
static Vec3 Boxminus(const Quat& q, const Quat& p){
  return Log(q*p.inverse());
}
static Mat3 SSM(const Vec3& vec){
  Mat3 mat;
  mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
  return mat;
}
static Mat3 RotMat(const Vec3& vec){
  const double a = vec.norm();
  if(a < 1e-8){
    return Mat3::Identity() + SSM(vec);
  } else{
    const Vec3 axis = vec.normalized();
    return Mat3::Identity() + sin(a)*SSM(axis) + (1-cos(a))*SSM(axis)*SSM(axis);
  }
}
static Mat3 GammaMat(const Vec3& vec){
  const double a = vec.norm();
  if(a < 1e-8){
    return Mat3::Identity() + 0.5 * SSM(vec);
  } else{
    const Vec3 axis = vec.normalized();
    return Mat3::Identity() + (1-cos(a))/a*SSM(axis) + (a-sin(a))/a*SSM(axis)*SSM(axis);
  }
}
static Mat3 FromTwoVectorsJac(const Vec3& a, const Vec3& b){
  const Vec3 cross = a.cross(b);
  const double crossNorm = cross.norm();
  Vec3 crossNormalized = cross/crossNorm;
  Mat3 crossNormalizedSqew = SSM(crossNormalized);
  const double c = a.dot(b);
  const double angle = std::acos(c);
  if(crossNorm<1e-6){
    if(c>0){
      return -SSM(b);
    } else {
      TSIF_LOGW("Warning: instable FromTwoVectorsJac!");
      return Mat3::Zero();
    }
  } else {
    return -1/crossNorm*(crossNormalized*b.transpose()-(crossNormalizedSqew*crossNormalizedSqew*SSM(b)*angle));
  }
}

typedef std::chrono::high_resolution_clock Clock;
typedef Clock::time_point TimePoint;
typedef Clock::duration Duration;
inline Duration fromSec(const double sec){
  return std::chrono::duration_cast < Duration
      > (std::chrono::duration<double>(sec));
}
inline double toSec(const Duration& duration){
  return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}
static std::string Print(TimePoint t){
  std::ostringstream out;
  out.precision(15);
  out << ((double)t.time_since_epoch().count()*Duration::period::num)/(Duration::period::den);
  return out.str();
}

/*! \brief Normal Random Number Generator
 *         Singleton class for generating normal random numbers (N(0,1)). Allows setting of seed.
 */
class NormalRandomNumberGenerator{
 public:
  void SetSeed(int s){
    generator_.seed(s);
  }
  double Get(){
    return distribution_(generator_);
  }
  template<int N>
  Vec<N> GetVec(){
    Vec<N> n;
    for(int i=0;i<N;i++){
      n(i) = Get();
    }
    return n;
  }
  static NormalRandomNumberGenerator& Instance(){
    static NormalRandomNumberGenerator instance;
    return instance;
  }
  std::default_random_engine& GetGenerator(){
    return generator_;
  }
 protected:
  std::default_random_engine generator_;
  std::normal_distribution<double> distribution_;
  NormalRandomNumberGenerator(): generator_(0), distribution_(0.0,1.0){
  }
};

template<typename T>
struct OptionLoaderTraits{
  static bool Get(T& x, const std::vector<std::string>& data){
    return false;
  }
};
template<>
struct OptionLoaderTraits<int>{
  static bool Get(int& x, const std::vector<std::string>& data){
    assert(data.size() == 1);
    x = stoi(data[0]);
    return true;
  }
};
template<>
struct OptionLoaderTraits<float>{
  static bool Get(float& x, const std::vector<std::string>& data){
    assert(data.size() == 1);
    x = stof(data[0]);
    return true;
  }
};
template<>
struct OptionLoaderTraits<double>{
  static bool Get(double& x, const std::vector<std::string>& data){
    assert(data.size() == 1);
    x = stod(data[0]);
    return true;
  }
};
template<int N>
struct OptionLoaderTraits<Vec<N>>{
  static bool Get(Vec<N>& x, const std::vector<std::string>& data){
    assert(data.size() == N);
    for(int i=0;i<N;i++){
      x(i) = stod(data[i]);
    }
    return true;
  }
};
template<>
struct OptionLoaderTraits<Quat>{
  static bool Get(Quat& x, const std::vector<std::string>& data){
    assert(data.size() == 4);
    x.w() = stod(data[0]);
    x.x() = stod(data[1]);
    x.y() = stod(data[2]);
    x.z() = stod(data[3]);
    return true;
  }
};

/*! \brief Option Loader
 *         Singleton class for interfacing option files.
 */
class OptionLoader{
 public:
  typedef std::vector<std::string> optionData;
  typedef std::map<std::string,optionData> FileData;
  std::map<std::string,FileData> data_;
  void LoadFile(const std::string& filename){
    if(data_.count(filename) == 0){
      FileData fileData;
      std::ifstream data(filename);
      std::string str;
      while(getline(data, str)){
        size_t prev = 0;
        std::vector<std::string> dataVector;
        bool isFirst = true;
        std::string name;
        while(prev <= str.size()){
          if(str[prev] == '#'){
            break;
          }
          const size_t next = str.find_first_of("\t ",prev);
          if(next>prev){
            if(isFirst){
              name = str.substr(prev,next-prev);
            } else {
              dataVector.push_back(str.substr(prev,next-prev));
            }
            isFirst = false;
          }
          prev = next != std::string::npos ? next + 1 : next;
        }
        if(!isFirst){
          fileData[name] = dataVector;
        }
      }
      data_[filename] = fileData;
    }
  }
  void PrintData(const FileData& d){
    for(const auto& entry : d){
      std::cout << entry.first << std::endl;
      for(const auto& value : entry.second){
        std::cout << value << "|";
      }
      std::cout << std::endl;
    }
  }
  template<typename T>
  bool Get(const std::string& filename, const std::string& name, T& x){
    LoadFile(filename);
    return OptionLoaderTraits<T>::Get(x,data_.at(filename).at(name));
  }
  template<typename T>
  T Get(const std::string& filename, const std::string& name){
    T x;
    Get(filename,name,x);
    return x;
  }
  static OptionLoader& Instance(){
    static OptionLoader instance;
    return instance;
  }
 protected:
  OptionLoader(){
  }
};

class Timer{
 public:
  Timer(){
    start_ = Clock::now();
    last_ = start_;
  }
  double GetIncr(){
    TimePoint now = Clock::now();
    double incr = toSec(now-last_);
    last_ = now;
    return incr;
  }
  double GetFull(){
    last_ = Clock::now();
    return toSec(last_-start_);
  }
  TimePoint start_;
  TimePoint last_;
};

} // namespace tsif

#endif /* TSIF_COMMON_HPP_ */
