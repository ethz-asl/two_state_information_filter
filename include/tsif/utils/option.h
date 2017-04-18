#ifndef TSIF_OPTION_HPP_
#define TSIF_OPTION_HPP_

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "tsif/utils/typedefs.h"

namespace tsif{

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

} // namespace tsif

#endif /* TSIF_OPTION_HPP_ */
