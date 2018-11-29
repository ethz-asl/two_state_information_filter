#ifndef TSIF_ELEMENT_VECTOR_H_
#define TSIF_ELEMENT_VECTOR_H_

#include "tsif/utils/common.h"
#include "tsif/element.h"

namespace tsif{

template <typename... Elements>
struct DimensionTrait;
template <>
struct DimensionTrait<>{
  static const int kDim = 0;
};
template <typename Element, typename... Elements>
struct DimensionTrait<Element,Elements...>{
  static const int kDim = Element::kDim + DimensionTrait<Elements...>::kDim;
};

template<bool... Bs>
struct CheckAllTrue;
template <bool B>
struct CheckAllTrue<B>{
  static const bool kVal = B;
};
template <bool B,bool... Bs>
struct CheckAllTrue<B,Bs...>{
  static const bool kVal = B & CheckAllTrue<Bs...>::kVal;
};


template<typename Derived, typename... Elements>
class ElementVectorBase{
 private:
  typedef std::tuple<Elements...> Tuple;

 public:
  static const int kN = sizeof...(Elements);
  static const bool kIsVectorSpace = CheckAllTrue<Elements::kIsVectorSpace...>::kVal;

  template<int I, int C = 0, typename std::enable_if<(std::tuple_element<C,Tuple>::type::kI == I)>::type* = nullptr>
  static constexpr int GetC(){
    return C;
  }
  template<int I, int C = 0, typename std::enable_if<((std::tuple_element<C,Tuple>::type::kI != I))>::type* = nullptr>
  static constexpr int GetC(){
    static_assert(C < kN-1, "Index not found in ElementVector");
    return GetC<I,C+1>();
  }

  template<int I>
  typename std::tuple_element<GetC<I>(),Tuple>::type::Type& Get(){
    return static_cast<Derived&>(*this).template Get<I>();
  }
  template<int I>
  const typename std::tuple_element<GetC<I>(),Tuple>::type::Type& Get() const{
    return static_cast<const Derived&>(*this).template Get<I>();
  }

  template<int I>
  typename std::tuple_element<GetC<I>(),Tuple>::type& GetElement(){
    return static_cast<Derived&>(*this).template GetElement<I>();
  }
  template<int I>
  const typename std::tuple_element<GetC<I>(),Tuple>::type& GetElement() const{
    return static_cast<const Derived&>(*this).template GetElement<I>();
  }

  template<typename Element, int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  static constexpr bool HasElement(){
    return std::is_same<typename std::tuple_element<C,Tuple>::type,Element>::value ? true : HasElement<Element,C+1>();
  }
  template<typename Element, int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  static constexpr bool HasElement(){
    return false;
  }

  template<int I, int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  static constexpr bool HasId(){
    return std::tuple_element<C,Tuple>::type::kI == I ? true : HasId<I,C+1>();
  }
  template<int I, int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  static constexpr bool HasId(){
    return false;
  }

  template<int C>
  static constexpr int GetId(){
    return std::tuple_element<C,Tuple>::type::kI;
  }

  template<int I>
  static constexpr int GetElementDim(){
    return std::tuple_element<GetC<I>(),Tuple>::type::kDim;
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  std::string Print() const{
    return GetElement<GetId<C>()>().Print() + "\n" + Print<C+1>();
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  std::string Print() const{
    return "";
  }

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void SetIdentity(){
    GetElement<GetId<C>()>().SetIdentity();
    SetIdentity<C+1>();
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void SetIdentity(){}

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void SetRandom(){
    GetElement<GetId<C>()>().SetRandom();
    SetRandom<C+1>();
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void SetRandom(){}

  template<typename OtherDerived, int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void Boxplus(const VecCRefX& vec, ElementVectorBase<OtherDerived,Elements...>& out) const{
    assert(vec.size() == Dim());
    typedef typename std::tuple_element<C,Tuple>::type E;
    if(E::kDim > 0){
      GetElement<E::kI>().Boxplus(vec.template block<E::kDim,1>(Start(E::kI),0),out.template GetElement<E::kI>());
    }
    Boxplus<OtherDerived,C+1>(vec,out);
  }
  template<typename OtherDerived, int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void Boxplus(const VecCRefX& vec, ElementVectorBase<OtherDerived,Elements...>& out) const{}

  template<typename OtherDerived, int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void Boxminus(const ElementVectorBase<OtherDerived,Elements...>& ref, VecRefX out) const{
    assert(out.size() == Dim());
    typedef typename std::tuple_element<C,Tuple>::type E;
    if(E::kDim > 0){
      GetElement<E::kI>().Boxminus(ref.template GetElement<E::kI>(),out.template block<E::kDim,1>(Start(E::kI),0));
    }
    Boxminus<OtherDerived,C+1>(ref,out);
  }
  template<typename OtherDerived, int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void Boxminus(const ElementVectorBase<OtherDerived,Elements...>& ref, VecRefX out) const{}

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void GetVec(VecRefX vec) const{
    assert(vec.size() == Dim());
    typedef typename std::tuple_element<C,Tuple>::type E;
    if(E::kDim>0){
      vec.template block<E::kDim,1>(Start(E::kI),0) = GetElement<E::kI>().GetVec();
    }
    GetVec<C+1>(vec);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void GetVec(VecRefX vec) const{}

  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  void Scale(double w){
    typedef typename std::tuple_element<C,Tuple>::type E;
    GetElement<E::kI>().Scale(w);
    Scale<C+1>(w);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  void Scale(double w){}

  int Start(int I) const{
    return static_cast<const Derived&>(*this).Start(I);
  }
  int Dim() const{
    return static_cast<const Derived&>(*this).Dim();
  }
};

template<typename... Elements>
class ElementVectorRef: public ElementVectorBase<ElementVectorRef<Elements...>,Elements...>{
 private:
  typedef ElementVectorBase<ElementVectorRef<Elements...>,Elements...> Base;
  typedef std::tuple<Elements...> Tuple;
  std::tuple<Elements&...> elements_;
  const std::map<int,int> startMap_;
  const int Dim_;

 public:
  static constexpr int kN = Base::kN;
  template<typename OtherDerived,typename... OtherElements>
  ElementVectorRef(ElementVectorBase<OtherDerived,OtherElements...>& elementVector):
    elements_(elementVector.template GetElement<Elements::kI>()...),
    startMap_{std::make_pair(Elements::kI,elementVector.Start(Elements::kI))...},
    Dim_(elementVector.Dim()){
  }
  ~ElementVectorRef(){}

  template<int I>
  typename std::tuple_element<Base::template GetC<I>(),Tuple>::type::Type& Get(){
    return std::get<this->template GetC<I>()>(elements_).Get();
  }
  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type::Type& Get() const{
    return std::get<this->template GetC<I>()>(elements_).Get();
  }
  template<int I>
  typename std::tuple_element<Base::template GetC<I>(),Tuple>::type& GetElement(){
    return std::get<this->template GetC<I>()>(elements_);
  }
  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type& GetElement() const{
    return std::get<this->template GetC<I>()>(elements_);
  }
  int Start(int I) const{
    return startMap_.at(I);
  }
  int Dim() const{
    return Dim_;
  }
};

template<typename... Elements>
class ElementVectorConstRef: public ElementVectorBase<ElementVectorConstRef<Elements...>,Elements...>{
 private:
  typedef ElementVectorBase<ElementVectorConstRef<Elements...>,Elements...> Base;
  typedef std::tuple<Elements...> Tuple;
  std::tuple<const Elements&...> elements_;
  const std::map<int,int> startMap_;
  const int Dim_;

 public:
  static constexpr int kN = Base::kN;
  template<typename OtherDerived,typename... OtherElements>
  ElementVectorConstRef(const ElementVectorBase<OtherDerived,OtherElements...>& elementVector):
    elements_(elementVector.template GetElement<Elements::kI>()...),
    startMap_{std::make_pair(Elements::kI,elementVector.Start(Elements::kI))...},
    Dim_(elementVector.Dim()){
  }
  ~ElementVectorConstRef(){}

  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type::Type& Get() const{
    return std::get<this->template GetC<I>()>(elements_).Get();
  }
  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type& GetElement() const{
    return std::get<this->template GetC<I>()>(elements_);
  }
  int Start(int I) const{
    return startMap_.at(I);
  }
  int Dim() const{
    return Dim_;
  }
};

template<typename... Elements>
class ElementVector: public ElementVectorBase<ElementVector<Elements...>,Elements...>{
 private:
  typedef ElementVectorBase<ElementVector<Elements...>,Elements...> Base;
  typedef std::tuple<Elements...> Tuple;
  alignas(16) Tuple elements_;

 public:
  typedef ElementVectorRef<Elements...> Ref;
  typedef ElementVectorConstRef<Elements...> CRef;
  static constexpr int kN = Base::kN;
  ElementVector(){}
  template<typename OtherDerived,typename... OtherElements>
  ElementVector(const ElementVectorBase<OtherDerived,OtherElements...>& elementVector):
    elements_(elementVector.template GetElement<Elements::kI>()...){}
  template<int C = 0, typename std::enable_if<(kN > C)>::type* = nullptr>
  ElementVector(const typename Elements::Type&... elements):
    elements_(elements...){}
  ~ElementVector(){}

  template<int I>
  typename std::tuple_element<Base::template GetC<I>(),Tuple>::type::Type& Get(){
    return std::get<this->template GetC<I>()>(elements_).Get();
  }
  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type::Type& Get() const{
    return std::get<this->template GetC<I>()>(elements_).Get();
  }

  template<int I>
  typename std::tuple_element<Base::template GetC<I>(),Tuple>::type& GetElement(){
    return std::get<this->template GetC<I>()>(elements_);
  }
  template<int I>
  const typename std::tuple_element<Base::template GetC<I>(),Tuple>::type& GetElement() const{
    return std::get<this->template GetC<I>()>(elements_);
  }
  static constexpr int Start(int I){
    return _Start<0>(I);
  }
  template<int C = 0, typename std::enable_if<(C < kN)>::type* = nullptr>
  static constexpr int _Start(int I){
    return (std::tuple_element<C,Tuple>::type::kI == I) ? 0 : std::tuple_element<C,Tuple>::type::kDim + _Start<C+1>(I);
  }
  template<int C = 0, typename std::enable_if<(C >= kN)>::type* = nullptr>
  static constexpr int _Start(int I){
    return -1;
  }
  static constexpr int Dim(){
    return DimensionTrait<Elements...>::kDim;
  }
};

template<typename ElementVector, typename Element, bool DoAddition = true>
struct AddElementTrait;
template <typename Element, typename... Elements>
struct AddElementTrait<ElementVector<Elements...>,Element,true>{
  typedef ElementVector<Elements...,Element> Type;
};
template <typename Element, typename... Elements>
struct AddElementTrait<ElementVector<Elements...>,Element,false>{
  typedef ElementVector<Elements...> Type;
};

template<typename ElementVectorA,typename ElementVectorB>
struct MergeTwo;

template<typename... Elements>
struct MergeTwo<ElementVector<Elements...>,ElementVector<>>{
  typedef ElementVector<Elements...> Type;
};
template<typename Element,typename... ElementsA,typename... ElementsB>
struct MergeTwo<ElementVector<ElementsA...>,ElementVector<Element,ElementsB...>>{
  static_assert(!ElementVector<ElementsA...>::template HasId<Element::kI>() |
                ElementVector<ElementsA...>::template HasElement<Element>(),"Conflicting Merge");
  typedef typename MergeTwo<typename AddElementTrait<ElementVector<ElementsA...>,Element,
                           !ElementVector<ElementsA...>::template HasElement<Element>()>::Type,
                            ElementVector<ElementsB...>>::Type Type;
};

template<typename... ElementVectors>
struct MergeTrait;
template<typename ElementVector>
struct MergeTrait<ElementVector>{
  typedef ElementVector Type;
};
template<typename ElementVectorA, typename ElementVectorB, typename... ElementVectors>
struct MergeTrait<ElementVectorA,ElementVectorB,ElementVectors...>{
  typedef typename MergeTrait<typename MergeTwo<ElementVectorA,ElementVectorB>::Type,
                              ElementVectors...>::Type Type;
};

class MeasEmpty: public ElementVector<>{
};

} // namespace tsif

#endif  // TSIF_ELEMENT_VECTOR_H_
