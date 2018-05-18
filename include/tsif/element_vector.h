#ifndef TSIF_ELEMENT_VECTOR_H_
#define TSIF_ELEMENT_VECTOR_H_

#include "tsif/utils/common.h"
#include "tsif/element.h"

namespace tsif{

template<int I, int... Range, int... Elements>
static int constexpr get_element_impl(std::integer_sequence<int, Range...>, std::integer_sequence<int, Elements...>){
  return (((Elements == I) * (Range + 1)) + ... + -1);
}

template<int I, int... Elements>
static int constexpr get_element(){
  return get_element_impl<I>(std::make_integer_sequence<int, sizeof...(Elements)>{}, std::integer_sequence<int, Elements...>{});
}

template<typename Derived, typename... Elements>
class ElementVectorBase{
 private:
  typedef std::tuple<Elements...> Tuple;

 public:
  static const int kN = sizeof...(Elements);
  static const bool kIsVectorSpace = (... & Elements::kIsVectorSpace);

  template<int I>
  static constexpr int GetC(){
    return get_element<I, Elements::kI...>();
  }

  template<int I>
  auto& Get(){
    return static_cast<Derived&>(*this).Get<I>();
  }
  template<int I>
  const auto& Get() const{
    return static_cast<const Derived&>(*this).Get<I>();
  }

  template<int I>
  auto& GetElement(){
    return static_cast<Derived&>(*this).GetElement<I>();
  }
  template<int I>
  const auto& GetElement() const{
    return static_cast<const Derived&>(*this).GetElement<I>();
  }

  template<typename Element>
  static constexpr bool HasElement(){
    return (std::is_same<Elements,Element>::value | ... | false);
  }

  template<int I>
  static constexpr bool HasId(){
    return ((Elements::kI == I) | ... | false);
  }

  template<int C>
  static constexpr int GetId(){
    return std::tuple_element<C,Tuple>::type::kI;
  }

  template<int I>
  static constexpr int GetElementDim(){
    return std::tuple_element<GetC<I>(),Tuple>::type::kDim;
  }

  std::string Print() const{
    return ((GetElement<Elements::kI>().Print() + "\n") + ... + "");
  }

  void SetIdentity(){
    (GetElement<Elements::kI>().SetIdentity(), ...);
  }

  void SetRandom(){
    (GetElement<Elements::kI>().SetRandom(), ...);
  }

  template<typename OtherDerived>
  void Boxplus(const VecCRefX& vec, ElementVectorBase<OtherDerived,Elements...>& out) const{
    assert(vec.size() == Dim());
    (GetElement<Elements::kI>().Boxplus(vec.template block<Elements::kDim,1>(Start(Elements::kI),0),out.GetElement<Elements::kI>()), ...);
  }

  template<typename OtherDerived>
  void Boxminus(const ElementVectorBase<OtherDerived,Elements...>& ref, VecRefX out) const{
    assert(out.size() == Dim());
    (GetElement<Elements::kI>().Boxminus(ref.GetElement<Elements::kI>(),out.template block<Elements::kDim,1>(Start(Elements::kI),0)), ...);
  }

  void GetVec(VecRefX vec) const{
    assert(vec.size() == Dim());
    (GetElement<Elements::kI>().GetVec(vec.template block<Elements::kDim,1>(Start(Elements::kI),0)), ...);
  }

  void Scale(double w){
    (GetElement<Elements::kI>().Scale(w), ...);
  }

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
    elements_(elementVector.GetElement<Elements::kI>()...),
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
    elements_(elementVector.GetElement<Elements::kI>()...),
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
    elements_(elementVector.GetElement<Elements::kI>()...){}
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
    return _Start(I);
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
    return (Elements::kDim + ... + 0);
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
