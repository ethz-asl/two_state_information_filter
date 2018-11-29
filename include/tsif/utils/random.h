#ifndef TSIF_RANDOM_HPP_
#define TSIF_RANDOM_HPP_

#include "tsif/utils/typedefs.h"

#include <random>

namespace tsif{

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

} // namespace tsif

#endif /* TSIF_RANDOM_HPP_ */
