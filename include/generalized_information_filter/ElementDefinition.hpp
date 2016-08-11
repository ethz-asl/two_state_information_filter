/*
 * ElementDefinition.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_ELEMENTDEFINITION_HPP_
#define GIF_ELEMENTDEFINITION_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/Element.hpp"

namespace GIF{

template<typename ElementType>
class ElementTraits{
 public:
  static constexpr int d_ = 0;
  static void print(const ElementType& x){}
  static void init(ElementType& x){}
  static void boxplus(const ElementType& in, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementType& out){
    out = in;
  }
  static void boxminus(const ElementType& in, const ElementType& ref, Eigen::Ref<Eigen::VectorXd> vec){}
};

class ElementDefinitionBase{
 public:
  ElementDefinitionBase(int d): d_(d){};
  virtual ~ElementDefinitionBase(){};
  inline int getDim(){
    return d_;
  }
  virtual ElementBase* newElement() = 0;
  virtual void print(const ElementBase* e) const = 0;
  virtual void init(ElementBase* e) = 0;
  virtual void boxplus(const ElementBase* ref, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementBase* out) = 0;
  virtual void boxminus(const ElementBase* in, const ElementBase* ref, Eigen::Ref<Eigen::VectorXd> vec) = 0;

 protected:
  const int d_;
};

template<typename ElementType>
class ElementDefinition: public ElementDefinitionBase{
 public:
  ElementDefinition(): ElementDefinitionBase(ElementTraits<ElementType>::d_){};
  Element<ElementType>* newElement(){
    return new Element<ElementType>(this);
  }
  void print(const ElementBase* e) const{
    ElementTraits<ElementType>::print((dynamic_cast<const Element<ElementType>*>(e))->get());
  }
  void init(ElementBase* e){
    ElementTraits<ElementType>::init((dynamic_cast<Element<ElementType>*>(e))->get());
  }
  void boxplus(const ElementBase* ref, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementBase* out){
    ElementTraits<ElementType>::boxplus(dynamic_cast<const Element<ElementType>*>(ref)->get(),vec,dynamic_cast<Element<ElementType>*>(out)->get());
  }
  void boxminus(const ElementBase* in, const ElementBase* ref, Eigen::Ref<Eigen::VectorXd> vec){
    ElementTraits<ElementType>::boxminus(dynamic_cast<const Element<ElementType>*>(in)->get(),dynamic_cast<const Element<ElementType>*>(ref)->get(),vec);
  }
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

#endif /* GIF_ELEMENTDEFINITION_HPP_ */
