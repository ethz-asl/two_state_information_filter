/*
 * ElementDefinition.hpp
 *
 *  Created on: Jul 29, 2016
 *      Author: Bloeschm
 */

#ifndef GIF_ELEMENTDEFINITION_HPP_
#define GIF_ELEMENTDEFINITION_HPP_

#include "generalized_information_filter/common.hpp"
#include "generalized_information_filter/State.hpp"

namespace GIF{

template<typename ElementType>
class PrintTraits{
 public:
  static void print(const ElementType& x){}
};

template<typename ElementType>
class InitTraits{
 public:
  static void init(ElementType& x){}
};

template<typename ElementType>
class DimTraits{
 public:
  static constexpr int d_ = 0;
};

template<typename ElementType>
class BoxplusTraits{
 public:
  static void boxplus(const ElementType& in, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementType& out){
    out = in;
  }
};

template<typename ElementType>
class BoxminusTraits{
 public:
  static void boxminus(const ElementType& in, const ElementType& ref, Eigen::Ref<Eigen::VectorXd> vec){}
};

class ElementDefinitionBase{
 public:
  ElementDefinitionBase(int d): d_(d){
    i_ = -1;
    j_ = -1;
  };
  ElementDefinitionBase(int d, int i, int j): d_(d), i_(i), j_(j){};
  virtual ~ElementDefinitionBase(){};
  inline int getIndex(){
    return j_;
  }
  inline int getDim(){
    return d_;
  }
  virtual ElementBase* newElement() = 0;
  virtual void print(const ElementBase* e) const = 0;
  virtual void init(ElementBase* e) = 0;
  virtual void boxplus(const ElementBase* ref, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementBase* out) = 0;
  virtual void boxminus(const ElementBase* in, const ElementBase* ref, Eigen::Ref<Eigen::VectorXd> vec) = 0;

 protected:
  int i_;
  int j_;
  const int d_;
};

template<typename ElementType>
class ElementDefinition: public ElementDefinitionBase{
 public:
  ElementDefinition(): ElementDefinitionBase(DimTraits<ElementType>::d_){};
  ElementDefinition(int i, int j): ElementDefinitionBase(DimTraits<ElementType>::d_,i,j){};
  ElementType& get(State* s){
    return (dynamic_cast<Element<ElementType>*>(s->getElement(i_)))->x_;
  }
  const ElementType& get(const State* s){
    return (dynamic_cast<const Element<ElementType>*>(s->getElement(i_)))->x_;
  }
  Element<ElementType>* newElement(){
    return new Element<ElementType>(this);
  }
  void print(const ElementBase* e) const{
    PrintTraits<ElementType>::print((dynamic_cast<const Element<ElementType>*>(e))->x_);
  }
  void init(ElementBase* e){
    InitTraits<ElementType>::init((dynamic_cast<Element<ElementType>*>(e))->x_);
  }
  void boxplus(const ElementBase* ref, const Eigen::Ref<const Eigen::VectorXd>& vec, ElementBase* out){
    BoxplusTraits<ElementType>::boxplus(dynamic_cast<const Element<ElementType>*>(ref)->x_,vec,dynamic_cast<Element<ElementType>*>(out)->x_);
  }
  void boxminus(const ElementBase* in, const ElementBase* ref, Eigen::Ref<Eigen::VectorXd> vec){
    BoxminusTraits<ElementType>::boxminus(dynamic_cast<const Element<ElementType>*>(in)->x_,dynamic_cast<const Element<ElementType>*>(ref)->x_,vec);
  }
};

// ==================== Traits Implementation ====================
template<>
class PrintTraits<double>{
 public:
  static void print(const double& x){
    std::cout << x << std::endl;
  }
};
template<int N>
class PrintTraits<Eigen::Matrix<double,N,1>>{
 public:
  static void print(const Eigen::Matrix<double,N,1>& x){
    std::cout << x.transpose() << std::endl;
  }
};
template<typename T, size_t N>
class PrintTraits<std::array<T,N>>{
 public:
  static void print(const std::array<T,N>& x){
    for(const T& i : x){
      PrintTraits<T>::print(i);
    }
  }
};

template<>
class InitTraits<double>{
 public:
  static void init(double& x){
    x = 0;
  }
};
template<int N>
class InitTraits<Eigen::Matrix<double,N,1>>{
 public:
  static void init(Eigen::Matrix<double,N,1>& x){
    x.setZero();
  }
};
template<typename T, size_t N>
class InitTraits<std::array<T,N>>{
 public:
  static void init(std::array<T,N>& x){
    for(T& i : x){
      InitTraits<T>::init(i);
    }
  }
};

template<>
class DimTraits<double>{
 public:
  static constexpr int d_ = 1;
};
template<int N>
class DimTraits<Eigen::Matrix<double,N,1>>{
 public:
  static constexpr int d_ = N;
};
template<typename T, size_t N>
class DimTraits<std::array<T,N>>{
 public:
  static constexpr int d_ = N*DimTraits<T>::d_;
};

template<>
class BoxplusTraits<double>{
 public:
  static void boxplus(const double& in, const Eigen::Ref<const Eigen::VectorXd>& vec, double& out){
    out = in+vec(0);
  }
};
template<int N>
class BoxplusTraits<Eigen::Matrix<double,N,1>>{
 public:
  static void boxplus(const Eigen::Matrix<double,N,1>& in, const Eigen::Ref<const Eigen::VectorXd>& vec, Eigen::Matrix<double,N,1>& out){
    out = in+vec;
  }
};
template<typename T, size_t N>
class BoxplusTraits<std::array<T,N>>{
 public:
  static void boxplus(const std::array<T,N>& in, const Eigen::Ref<const Eigen::VectorXd>& vec, std::array<T,N>& out){
    for(int i=0;i<N;i++){
      BoxplusTraits<T>::boxplus(in[i],vec.block<DimTraits<T>::d_,1>(i*DimTraits<T>::d_,0),out[i]);
    }
  }
};

template<>
class BoxminusTraits<double>{
 public:
  static void boxminus(const double& in, const double& ref, Eigen::Ref<Eigen::VectorXd> vec){
    vec(0) = in-ref;
  }
};
template<int N>
class BoxminusTraits<Eigen::Matrix<double,N,1>>{
 public:
  static void boxminus(const Eigen::Matrix<double,N,1>& in, const Eigen::Matrix<double,N,1>& ref, Eigen::Ref<Eigen::VectorXd> vec){
    vec = in-ref;
  }
};
template<typename T, size_t N>
class BoxminusTraits<std::array<T,N>>{
 public:
  static void boxminus(const std::array<T,N>& in, const std::array<T,N>& ref, Eigen::Ref<Eigen::VectorXd> vec){
    for(int i=0;i<N;i++){
      BoxminusTraits<T>::boxminus(in[i],ref[i],vec.block<DimTraits<T>::d_,1>(i*DimTraits<T>::d_,0));
    }
  }
};

}

#endif /* GIF_ELEMENTDEFINITION_HPP_ */
