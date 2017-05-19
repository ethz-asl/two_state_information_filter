#ifndef TSIF_ELEMENT_H_
#define TSIF_ELEMENT_H_

#include "tsif/utils/common.h"
#include "unit_vector.h"

namespace tsif{

template<typename T>
struct ElementTraits{
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = 0;
  static std::string Print(const T& x){
    return "";
  }
  static T Identity(){
    T x;
    return x;
  }
  static void SetRandom(T& x){
  }
  static void Boxplus(const T& ref, const VecCRef<kDim>& vec, T& out){
    out = ref;   // Must be computable in-place
  }
  static void Boxminus(const T& in, const T& ref, VecRef<kDim> vec){
  }
  static Mat<kDim, kDim> BoxplusJacInp(const T& in, const VecCRef<kDim>& vec){
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxplusJacVec(const T& in, const VecCRef<kDim>& vec){
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxminusJacInp(const T& in, const T& ref){
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxminusJacRef(const T& in, const T& ref){
    return Mat<kDim, kDim>::Identity();
  }
  static Vec<kDim> GetVec(const T& x){
    return Vec<kDim>::Zero();
  }
  static void Scale(double w, T& x){
  }
};

template<typename T, int I>
class Element{
 private:
	alignas(16) T x_;
  typedef ElementTraits<T> Traits;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef T Type;
  static const int kI = I;
  static const int kDim = Traits::kDim*(int)(I>=0);
  static const bool kIsVectorSpace = Traits::kIsVectorSpace;
  Element(){}
  Element(const T& x): x_(x){}
  T& Get(){
    return x_;
  }
  const T& Get() const{
    return x_;
  }
  std::string Print() const{
    std::ostringstream out;
    out << (I>=0 ? "\033[32m" : "\033[33m");
    out << Traits::Print(x_) << "\033[0m";
    return out.str();
  }
  void SetIdentity(){
    x_ = Traits::Identity();
  }
  void SetRandom(){
    Traits::SetRandom(x_);
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  void Boxplus(const VecCRef<kDim>& vec, Element<T,I>& out) const{
    Traits::Boxplus(x_, vec, out.Get());
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  void Boxplus(const VecCRef<kDim>& vec, Element<T,I>& out) const{
    out.Get() = x_;
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  void Boxminus(const Element<T,I>& ref, VecRef<kDim> vec) const{
    Traits::Boxminus(x_, ref.Get(), vec);
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  void Boxminus(const Element<T,I>& ref, VecRef<kDim> vec) const{
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  Mat<kDim> BoxplusJacInp(const VecCRef<kDim>& vec) const{
    return Traits::BoxplusJacInp(x_, vec);
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  Mat<kDim> BoxplusJacInp(const VecCRef<kDim>& vec) const{
    return Mat<kDim>::Zero();
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  Mat<kDim> BoxplusJacVec(const VecCRef<kDim>& vec) const{
    return Traits::BoxplusJacVec(x_, vec);
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  Mat<kDim> BoxplusJacVec(const VecCRef<kDim>& vec) const{
    return Mat<kDim>::Zero();
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  Mat<kDim> BoxminusJacInp(const Element<T,I>& ref) const{
    return Traits::BoxminusJacInp(x_, ref.Get());
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  Mat<kDim> BoxminusJacInp(const Element<T,I>& ref) const{
    return Mat<kDim>::Zero();
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  Mat<kDim> BoxminusJacRef(const Element<T,I>& ref) const{
    return Traits::BoxminusJacRef(x_, ref.Get());
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  Mat<kDim> BoxminusJacRef(const Element<T,I>& ref) const{
    return Mat<kDim>::Zero();
  }
  template<bool B = I>=0,typename std::enable_if<B>::type* = nullptr>
  Vec<kDim> GetVec() const{
    return Traits::GetVec(x_);
  }
  template<bool B = I>=0,typename std::enable_if<!B>::type* = nullptr>
  Vec<kDim> GetVec() const{
    return Vec<kDim>::Zero();
  }
  void Scale(double w){
    Traits::Scale(w,x_);
  }
};

template<typename T, int I>
const int Element<T,I>::kI;

// ==================== Traits Implementation ==================== //
/*! \brief Scalar Trait.
 *         Element trait for regular scalar.
 */
template<>
class ElementTraits<double>{
 public:
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = 1;
  static std::string Print(const double& x){
    std::ostringstream out;
    out << x;
    return out.str();
  }
  static double Identity(){
    return 0;
  }
  static void SetRandom(double& x){
    x = NormalRandomNumberGenerator::Instance().Get();
  }
  static void Boxplus(const double& in, const VecCRef<kDim>& vec, double& out){
    out = in + vec(0);
  }
  static void Boxminus(const double& in, const double& ref, VecRef<kDim> vec){
    vec(0) = in - ref;
  }
  static Mat<kDim> BoxplusJacInp(const double& in, const VecCRef<kDim>& vec){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxplusJacVec(const double& in, const VecCRef<kDim>& vec){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacInp(const double& in, const double& ref){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacRef(const double& in, const double& ref){
    return -Mat<kDim>::Identity();
  }
  static Vec<kDim> GetVec(const double& x){
    return Vec<1>(x);
  }
  static void Scale(double w, double& x){
    x *= w;
  }
};

/*! \brief Vector Trait.
 *         Element trait for vector of scalars
 */
template<int N>
class ElementTraits<Vec<N>>{
 public:
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = N;
  static std::string Print(const Vec<N>& x){
    std::ostringstream out;
    out << x.transpose();
    return out.str();
  }
  static Vec<N> Identity(){
    return Vec<N>::Zero();
  }
  static void SetRandom(Vec<N>& x){
    x = NormalRandomNumberGenerator::Instance().GetVec<N>();
  }
  static void Boxplus(const Vec<N>& in, const VecCRef<kDim>& vec, Vec<N>& out){
    out = in + vec;
  }
  static void Boxminus(const Vec<N>& in, const Vec<N>& ref, VecRef<kDim> vec){
    vec = in - ref;
  }
  static Mat<kDim> BoxplusJacInp(const Vec<N>& in, const VecCRef<kDim>& vec){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxplusJacVec(const Vec<N>& in, const VecCRef<kDim>& vec){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacInp(const Vec<N>& in, const Vec<N>& ref){
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacRef(const Vec<N>& in, const Vec<N>& ref){
    return -Mat<kDim>::Identity();
  }
  static Vec<kDim> GetVec(const Vec<N>& x){
    return x;
  }
  static void Scale(double w, Vec<N>& x){
    x *= w;
  }
};

/*! \brief Unit Quaternion Trait.
 *         Element trait for unit quaternion. Employed to represent
 *         orientations.
 */
template<>
class ElementTraits<Quat>{
 public:
  static constexpr bool kIsVectorSpace = false;
  static constexpr int kDim = 3;
  static std::string Print(const Quat& x){
    std::ostringstream out;
    out << x.w() << " " << x.x() << " " << x.y() << " " << x.z();
    return out.str();
  }
  static Quat Identity(){
    return Quat::Identity();
  }
  static void SetRandom(Quat& x){
    x.w() = NormalRandomNumberGenerator::Instance().Get();
    x.x() = NormalRandomNumberGenerator::Instance().Get();
    x.y() = NormalRandomNumberGenerator::Instance().Get();
    x.z() = NormalRandomNumberGenerator::Instance().Get();
    x.normalize();
  }
  static void Boxplus(const Quat& in, const VecCRef<kDim>& vec, Quat& out){
    out = tsif::Boxplus(in,vec);
    out.normalize();
  }
  static void Boxminus(const Quat& in, const Quat& ref, VecRef<kDim> vec){
    vec = tsif::Boxminus(in,ref);
  }
  static Mat<kDim> BoxplusJacInp(const Quat& in, const VecCRef<kDim>& vec){
    return Exp(vec).toRotationMatrix();
  }
  static Mat<kDim> BoxplusJacVec(const Quat& in, const VecCRef<kDim>& vec){
    return GammaMat(vec);
  }
  static Mat<kDim> BoxminusJacInp(const Quat& in, const Quat& ref){
    return GammaMat(tsif::Boxminus(in,ref)).inverse();
  }
  static Mat<kDim> BoxminusJacRef(const Quat& in, const Quat& ref){
    return -GammaMat(tsif::Boxminus(in,ref)).inverse() * (in*ref.inverse()).toRotationMatrix();
  }
  static Vec<kDim> GetVec(const Quat& x){
    return Log(x);
  }
  static void Scale(double w, Quat& x){
    x = Exp(w*Log(x));
  }
};

/*! \brief Unit Vector Trait.
 *         Element trait for unit vectors.
 */
template<>
class ElementTraits<UnitVector>{
 public:
  static constexpr bool kIsVectorSpace = false;
  static constexpr int kDim = 2;
  static std::string Print(const UnitVector& x){
    std::ostringstream out;
    out << x.GetVec().transpose();
    return out.str();
  }
  static UnitVector Identity(){
    return UnitVector(Quat::Identity());
  }
  static void SetRandom(UnitVector& x){
    x.SetRandom();
  }
  static void Boxplus(const UnitVector& in, const VecCRef<kDim>& vec, UnitVector& out){
    in.Boxplus(vec,out);
  }
  static void Boxminus(const UnitVector& in, const UnitVector& ref, VecRef<kDim> vec){
    in.Boxminus(ref,vec);
  }
  static Mat<kDim> BoxplusJacInp(const UnitVector& in, const VecCRef<kDim>& vec){
    Mat<kDim> J;
    in.BoxplusJacInp(vec,J);
    return J;
  }
  static Mat<kDim> BoxplusJacVec(const UnitVector& in, const VecCRef<kDim>& vec){
    Mat<kDim> J;
    in.BoxplusJacVec(vec,J);
    return J;
  }
  static Mat<kDim> BoxminusJacInp(const UnitVector& in, const UnitVector& ref){
    Mat<kDim> J;
    in.BoxminusJacInp(ref,J);
    return J;
  }
  static Mat<kDim> BoxminusJacRef(const UnitVector& in, const UnitVector& ref){
    Mat<kDim> J;
    in.BoxminusJacRef(ref,J);
    return J;
  }
  static Vec<kDim> GetVec(const UnitVector& x){
    Vec<kDim> vec;
    x.Boxminus(Identity(),vec);
    return vec;
  }
  static void Scale(double w, UnitVector& x){
    Vec<kDim> vec;
    x.Boxminus(Identity(),vec);
    vec *= w;
    Identity().Boxplus(vec,x);
  }
};

/*! \brief Array Trait.
 *         Element trait for array of specific sub-elements
 */
template<typename T, size_t N>
class ElementTraits<std::array<T, N>>{
 public:
  typedef ElementTraits<T> Traits;
  using array = std::array<T, N>;
  static constexpr bool kIsVectorSpace = Traits::kIsVectorSpace;
  static constexpr int kElementDim = Traits::kDim;
  static constexpr int kDim = N * kElementDim;
  static std::string Print(const array& x){
    std::ostringstream out;
    for (const T& i : x){
      out << Traits::Print(i) << "\t";
    }
    return out.str();

  }
  static array Identity(){
    array x;
    for (T& i : x){
      i = Traits::Identity();
    }
    return x;
  }
  static void SetRandom(array& x){
    for (T& i : x){
      Traits::SetRandom(i);
    }
  }
  static void Boxplus(const array& in, const VecCRef<kDim>& vec, array& out){
    for (int i = 0; i < N; i++){
      Traits::Boxplus(in.at(i), vec.template block<kElementDim, 1>(i * kElementDim, 0), out.at(i));
    }
  }
  static void Boxminus(const array& in, const array& ref, VecRef<kDim> vec){
    for (int i = 0; i < N; i++){
      Traits::Boxminus(in.at(i), ref.at(i), vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
  }
  static Mat<kDim> BoxplusJacInp(const array& in, const VecCRef<kDim>& vec){
    Mat<kDim> J = Mat<kDim>::Zero();
    for (int i = 0; i < N; i++){
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxplusJacInp(in.at(i), vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
    return J;
  }
  static Mat<kDim> BoxplusJacVec(const array& in, const VecCRef<kDim>& vec){
    Mat<kDim> J = Mat<kDim>::Zero();
    for (int i = 0; i < N; i++){
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxplusJacVec(in.at(i), vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
    return J;
  }
  static Mat<kDim> BoxminusJacInp(const array& in, const array& ref){
    Mat<kDim> J = Mat<kDim>::Zero();
    for (int i = 0; i < N; i++){
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxminusJacInp(in.at(i), ref.at(i));
    }
    return J;
  }
  static Mat<kDim> BoxminusJacRef(const array& in, const array& ref){
    Mat<kDim> J = Mat<kDim>::Zero();
    for (int i = 0; i < N; i++){
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxminusJacRef(in.at(i), ref.at(i));
    }
    return J;
  }
  static Vec<kDim> GetVec(const array& x){
    Vec<kDim> vec;
    for (int i = 0; i < N; i++){
      vec.template segment<kElementDim>(i * kElementDim) = x.at(i);
    }
    return vec;
  }
  static void Scale(double w, array& x){
    for (int i = 0; i < N; i++){
      Traits::Scale(w,x.at(i));
    }
  }
};

} // namespace tsif

#endif  // TSIF_ELEMENT_H_
