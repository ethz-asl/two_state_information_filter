/*
 * Element.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_ELEMENT_HPP_
#define GIF_ELEMENT_HPP_

#include "generalized_information_filter/common.hpp"

namespace GIF{

template<typename T>
class ElementTraits{
 public:
  static constexpr int d_ = 0;
  static void print(const T& x){}
  static void init(T& x){}
  static void boxplus(const T& in, const Eigen::Ref<const Eigen::VectorXd>& vec, T& out){
    out = in;
  }
  static void boxminus(const T& in, const T& ref, Eigen::Ref<Eigen::VectorXd> vec){}
};

class ElementBase{
 public:
  ElementBase(){};
  virtual ~ElementBase(){};
  virtual ElementBase& operator=(const ElementBase& other) = 0;
  virtual int getDim() const = 0;
  virtual void print() const = 0;
  virtual void init() = 0;
  virtual void boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec, const std::shared_ptr<ElementBase>& out) const = 0;
  virtual void boxminus(const std::shared_ptr<const ElementBase>& ref, Eigen::Ref<Eigen::VectorXd> vec)  const= 0;
};

template<typename T>
class ElementDefinition;

template<typename T>
class Element: public ElementBase{
 public:
  Element(const ElementDefinition<T>* def): def_(def){};
  virtual ~Element(){};
  Element<T>& operator=(const Element<T>& other){
    get() = other.get();
    return *this;
  }
  virtual ElementBase& operator=(const ElementBase& other){
    *this = dynamic_cast<const Element<T>&>(other);
    return *this;
  }
  int getDim() const{
    return def_->getDim();
  }
  void print() const{
    ElementTraits<T>::print(get());
  }
  void init(){
    ElementTraits<T>::init(get());
  }
  void boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec, const std::shared_ptr<ElementBase>& out) const{
    ElementTraits<T>::boxplus(get(),vec,std::dynamic_pointer_cast<Element<T>>(out)->get());
  }
  void boxminus(const std::shared_ptr<const ElementBase>& ref, Eigen::Ref<Eigen::VectorXd> vec) const{
    ElementTraits<T>::boxminus(get(),std::dynamic_pointer_cast<const Element<T>>(ref)->get(),vec);
  }
  T& get(){
    return x_;
  }
  const T& get() const{
    return x_;
  }
 protected:
  T x_;
  const ElementDefinition<T>* def_;
};

// ==================== Traits Implementation ====================
template<>
class ElementTraits<double>{
 public:
  static constexpr int d_ = 1;
  static void print(const double& x){
    std::cout << x << std::endl;
  }
  static void init(double& x){
    x = 0;
  }
  static void boxplus(const double& in, const Eigen::Ref<const Eigen::VectorXd>& vec, double& out){
    out = in+vec(0);
  }
  static void boxminus(const double& in, const double& ref, Eigen::Ref<Eigen::VectorXd> vec){
    vec(0) = in-ref;
  }
};
template<int N>
class ElementTraits<Eigen::Matrix<double,N,1>>{
 public:
  static constexpr int d_ = N;
  static void print(const Eigen::Matrix<double,N,1>& x){
    std::cout << x.transpose() << std::endl;
  }
  static void init(Eigen::Matrix<double,N,1>& x){
    x.setZero();
  }
  static void boxplus(const Eigen::Matrix<double,N,1>& in, const Eigen::Ref<const Eigen::VectorXd>& vec, Eigen::Matrix<double,N,1>& out){
    out = in+vec;
  }
  static void boxminus(const Eigen::Matrix<double,N,1>& in, const Eigen::Matrix<double,N,1>& ref, Eigen::Ref<Eigen::VectorXd> vec){
    vec = in-ref;
  }
};
template<typename T, size_t N>
class ElementTraits<std::array<T,N>>{
 public:
  static constexpr int d_ = N*ElementTraits<T>::d_;
  static void print(const std::array<T,N>& x){
    for(const T& i : x){
      ElementTraits<T>::print(i);
    }
  }
  static void init(std::array<T,N>& x){
    for(T& i : x){
      ElementTraits<T>::init(i);
    }
  }
  static void boxplus(const std::array<T,N>& in, const Eigen::Ref<const Eigen::VectorXd>& vec, std::array<T,N>& out){
    for(int i=0;i<N;i++){
      ElementTraits<T>::boxplus(in[i],vec.block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0),out[i]);
    }
  }
  static void boxminus(const std::array<T,N>& in, const std::array<T,N>& ref, Eigen::Ref<Eigen::VectorXd> vec){
    for(int i=0;i<N;i++){
      ElementTraits<T>::boxminus(in[i],ref[i],vec.block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0));
    }
  }
};

}

#endif /* GIF_ELEMENT_HPP_ */
