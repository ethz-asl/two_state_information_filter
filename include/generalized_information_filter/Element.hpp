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
class ElementTraits{ // Default implementation for zero dimension elements (may hold data which is not actively estimated/optimized)
 public:
  static constexpr int d_ = 0;
  static void print(const T& x){}
  static const T identity(){
    T x;
    return x;
  }
  static void setIdentity(T& x){}
  static void setRandom(T& x, int& s){}
  static void boxplus(const T& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec, T& out){
    out = in;
  }
  static void boxminus(const T& in, const T& ref, Eigen::Ref<Eigen::Matrix<double,d_,1>> vec){} // Must be computable in-place
  static Eigen::Matrix<double,d_,d_> boxplusJacInp(const T& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacVec(const T& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacInp(const T& in, const T& ref){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacRef(const T& in, const T& ref){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
};

class ElementBase{
 public:
  ElementBase(){};
  virtual ~ElementBase(){};
  virtual ElementBase& operator=(const ElementBase& other) = 0;
  virtual int getDim() const = 0;
  virtual void print() const = 0;
  virtual void setIdentity() = 0;
  virtual void setRandom(int& s) = 0;
  virtual void boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec, const std::shared_ptr<ElementBase>& out) const = 0;
  virtual void boxminus(const std::shared_ptr<const ElementBase>& ref, Eigen::Ref<Eigen::VectorXd> vec) const = 0;
  virtual MXD boxplusJacInp(const Eigen::Ref<const Eigen::VectorXd>& vec) const = 0;
  virtual MXD boxplusJacVec(const Eigen::Ref<const Eigen::VectorXd>& vec) const = 0;
  virtual MXD boxminusJacInp(const std::shared_ptr<ElementBase>& ref) const = 0;
  virtual MXD boxminusJacRef(const std::shared_ptr<ElementBase>& ref) const = 0;
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
  inline int getDim() const{
    return def_->getDim();
  }
  void print() const{
    ElementTraits<T>::print(get());
  }
  void setIdentity(){
    ElementTraits<T>::setIdentity(get());
  }
  void setRandom(int& s){
    ElementTraits<T>::setRandom(get(),s);
  }
  void boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec, const std::shared_ptr<ElementBase>& out) const{
    ElementTraits<T>::boxplus(get(),vec,std::dynamic_pointer_cast<Element<T>>(out)->get());
  }
  void boxminus(const std::shared_ptr<const ElementBase>& ref, Eigen::Ref<Eigen::VectorXd> vec) const{
    ElementTraits<T>::boxminus(get(),std::dynamic_pointer_cast<const Element<T>>(ref)->get(),vec);
  }
  MXD boxplusJacInp(const Eigen::Ref<const Eigen::VectorXd>& vec) const{
    return ElementTraits<T>::boxplusJacInp(get(),vec);
  }
  MXD boxplusJacVec(const Eigen::Ref<const Eigen::VectorXd>& vec) const{
    return ElementTraits<T>::boxplusJacVec(get(),vec);
  }
  MXD boxminusJacInp(const std::shared_ptr<ElementBase>& ref) const{
    return ElementTraits<T>::boxminusJacInp(get(),std::dynamic_pointer_cast<const Element<T>>(ref)->get());
  }
  MXD boxminusJacRef(const std::shared_ptr<ElementBase>& ref) const{
    return ElementTraits<T>::boxminusJacRef(get(),std::dynamic_pointer_cast<const Element<T>>(ref)->get());
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


// ==================== Traits Implementation ==================== //
template<>
class ElementTraits<double>{
 public:
  static constexpr int d_ = 1;
  static void print(const double& x){
    std::cout << x << std::endl;
  }
  static const double identity(){
    return 0;
  }
  static void setIdentity(double& x){
    x = 0;
  }
  static void setRandom(double& x, int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    x = distribution(generator);
    ++s;
  }
  static void boxplus(const double& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec, double& out){
    out = in+vec(0);
  }
  static void boxminus(const double& in, const double& ref, Eigen::Ref<Eigen::Matrix<double,d_,1>> vec){
    vec(0) = in-ref;
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacInp(const double& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacVec(const double& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacInp(const double& in, const double& ref){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacRef(const double& in, const double& ref){
    return -Eigen::Matrix<double,d_,d_>::Identity();
  }
};
template<int N>
class ElementTraits<Eigen::Matrix<double,N,1>>{
 public:
  static constexpr int d_ = N;
  static void print(const Eigen::Matrix<double,N,1>& x){
    std::cout << x.transpose() << std::endl;
  }
  static const Eigen::Matrix<double,N,1> identity(){
    return Eigen::Matrix<double,N,1>::Zero();
  }
  static void setIdentity(Eigen::Matrix<double,N,1>& x){
    x.setZero();
  }
  static void setRandom(Eigen::Matrix<double,N,1>& x, int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    for(unsigned int i=0;i<N;i++){
      x(i) = distribution(generator);
    }
    ++s;
  }
  static void boxplus(const Eigen::Matrix<double,N,1>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec, Eigen::Matrix<double,N,1>& out){
    out = in+vec;
  }
  static void boxminus(const Eigen::Matrix<double,N,1>& in, const Eigen::Matrix<double,N,1>& ref, Eigen::Ref<Eigen::Matrix<double,d_,1>> vec){
    vec = in-ref;
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacInp(const Eigen::Matrix<double,N,1>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacVec(const Eigen::Matrix<double,N,1>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacInp(const Eigen::Matrix<double,N,1>& in, const Eigen::Matrix<double,N,1>& ref){
    return Eigen::Matrix<double,d_,d_>::Identity();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacRef(const Eigen::Matrix<double,N,1>& in, const Eigen::Matrix<double,N,1>& ref){
    return -Eigen::Matrix<double,d_,d_>::Identity();
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
  static const std::array<T,N> identity(){
    std::array<T,N> x;
    setIdentity(x);
    return x;
  }
  static void setIdentity(std::array<T,N>& x){
    for(T& i : x){
      ElementTraits<T>::setIdentity(i);
    }
  }
  static void setRandom(std::array<T,N>& x, int& s){
    for(T& i : x){
      ElementTraits<T>::setRandom(i, s);
    }
  }
  static void boxplus(const std::array<T,N>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec, std::array<T,N>& out){
    for(int i=0;i<N;i++){
      ElementTraits<T>::boxplus(in[i],vec.template block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0),out[i]);
    }
  }
  static void boxminus(const std::array<T,N>& in, const std::array<T,N>& ref, Eigen::Ref<Eigen::Matrix<double,d_,1>> vec){
    for(int i=0;i<N;i++){
      ElementTraits<T>::boxminus(in[i],ref[i],vec.template block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0));
    }
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacInp(const std::array<T,N>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    Eigen::Matrix<double,d_,d_> J;
    J.setZero();
    for(int i=0;i<N;i++){
      J.template block<ElementTraits<T>::d_,ElementTraits<T>::d_>(i*ElementTraits<T>::d_,i*ElementTraits<T>::d_) = ElementTraits<T>::boxplusJacInp(in[i],vec.template block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0));
    }
    return J;
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacVec(const std::array<T,N>& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    Eigen::Matrix<double,d_,d_> J;
    J.setZero();
    for(int i=0;i<N;i++){
      J.template block<ElementTraits<T>::d_,ElementTraits<T>::d_>(i*ElementTraits<T>::d_,i*ElementTraits<T>::d_) = ElementTraits<T>::boxplusJacVec(in[i],vec.template block<ElementTraits<T>::d_,1>(i*ElementTraits<T>::d_,0));
    }
    return J;
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacInp(const std::array<T,N>& in, const std::array<T,N>& ref){
    Eigen::Matrix<double,d_,d_> J;
    J.setZero();
    for(int i=0;i<N;i++){
      J.template block<ElementTraits<T>::d_,ElementTraits<T>::d_>(i*ElementTraits<T>::d_,i*ElementTraits<T>::d_) = ElementTraits<T>::boxminusJacInp(in[i],ref[i]);
    }
    return J;
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacRef(const std::array<T,N>& in, const std::array<T,N>& ref){
    Eigen::Matrix<double,d_,d_> J;
    J.setZero();
    for(int i=0;i<N;i++){
      J.template block<ElementTraits<T>::d_,ElementTraits<T>::d_>(i*ElementTraits<T>::d_,i*ElementTraits<T>::d_) = ElementTraits<T>::boxminusJacRef(in[i],ref[i]);
    }
    return J;
  }
};
template<>
class ElementTraits<QPD>{
 public:
  static constexpr int d_ = 3;
  static void print(const QPD& x){
    std::cout << x << std::endl;
  }
  static const QPD identity(){
    return QPD();
  }
  static void setIdentity(QPD& x){
    x.setIdentity();
  }
  static void setRandom(QPD& x, int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    x.toImplementation().w() = distribution(generator);
    x.toImplementation().x() = distribution(generator);
    x.toImplementation().y() = distribution(generator);
    x.toImplementation().z() = distribution(generator);
    x.fix();
    ++s;
  }
  static void boxplus(const QPD& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec, QPD& out){
    out = in.boxPlus(vec);
  }
  static void boxminus(const QPD& in, const QPD& ref, Eigen::Ref<Eigen::Matrix<double,d_,1>> vec){
    vec = in.boxMinus(ref);
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacInp(const QPD& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    MPD m = m.exponentialMap(vec);
    return m.matrix();
  }
  static Eigen::Matrix<double,d_,d_> boxplusJacVec(const QPD& in, const Eigen::Ref<const Eigen::Matrix<double,d_,1>>& vec){
    return Lmat(vec);
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacInp(const QPD& in, const QPD& ref){
    return Lmat(in.boxMinus(ref)).inverse();
  }
  static Eigen::Matrix<double,d_,d_> boxminusJacRef(const QPD& in, const QPD& ref){
    return -Lmat(in.boxMinus(ref)).inverse()*MPD(in*ref.inverted()).matrix();
  }
};

}

#endif /* GIF_ELEMENT_HPP_ */
