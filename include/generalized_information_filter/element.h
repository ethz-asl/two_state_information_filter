#ifndef GIF_ELEMENT_HPP_
#define GIF_ELEMENT_HPP_

#include "generalized_information_filter/common.h"

namespace GIF {

template<typename T>
class Element;
template<typename T>
class ElementDescription;

/*! \brief Element traits.
 *         Default implementation for zero dimension elements,
 *         may hold data which is not actively estimated/optimized.
 */
template<typename T>
struct ElementTraits {
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = 0;
  static void Print(const T& x) {
  }
  static const T Identity() {
    T x;
    return x;
  }
  static void SetIdentity(T& x) {
  }
  static void SetRandom(T& x, int& s) {
  }
  static void Boxplus(const T& ref, const VecCRef<kDim>& vec, T& out) {
    out = ref;   // Must be computable in-place
  }
  static void Boxminus(const T& in, const T& ref, VecRef<kDim> vec) {
  }
  static Mat<kDim, kDim> BoxplusJacInp(const T& in, const VecCRef<kDim>& vec) {
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxplusJacVec(const T& in, const VecCRef<kDim>& vec) {
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxminusJacInp(const T& in, const T& ref) {
    return Mat<kDim, kDim>::Identity();
  }
  static Mat<kDim, kDim> BoxminusJacRef(const T& in, const T& ref) {
    return Mat<kDim, kDim>::Identity();
  }
};

/*! \brief Element Base.
 *         Will be used to store a vector of elements in ElementVector
 */
class ElementBase {
 public:
  typedef std::shared_ptr<ElementBase> Ptr;
  typedef std::shared_ptr<const ElementBase> CPtr;
  ElementBase() {
  }
  virtual ~ElementBase() {
  }
  virtual ElementBase& operator=(const ElementBase& other) = 0;
  virtual int GetDim() const = 0;
  virtual void Print() const = 0;
  virtual void SetIdentity() = 0;
  virtual void SetRandom(int& s) = 0;
  virtual void Boxplus(const VecCRefX& vec,ElementBase* out) const = 0;
  virtual void Boxminus(const ElementBase& ref, VecRefX vec) const = 0;
  virtual MatX BoxplusJacInp(const VecCRefX& vec) const = 0;
  virtual MatX BoxplusJacVec(const VecCRefX& vec) const = 0;
  virtual MatX BoxminusJacInp(const ElementBase& ref) const = 0;
  virtual MatX BoxminusJacRef(const ElementBase& ref) const = 0;
  template<typename T>
  T& GetValue() {
    Element<T>* const element =  dynamic_cast<Element<T>*>(this);
    assert(element != nullptr);
    return element->GetValue();
  }
  template<typename T>
  const T& GetValue() const {
    const Element<T>* const element =  dynamic_cast<const Element<T>*>(this);
    assert(element != nullptr);
    return element->GetValue();
  }
};

/*! \brief Templated form of Element.
 *         Forwards the virtual methods of the base class to the corresponding
 *         trait implementation.
 */
template<typename T>
class Element : public ElementBase {
 public:
  typedef ElementTraits<T> Traits;
  Element(const ElementDescription<T>* description)
      : description_(description) {
    assert(description != nullptr);
  }
  virtual ~Element() {
  }
  virtual Element<T>& operator=(const Element<T>& other) {
    GetValue() = other.GetValue();
    return *this;
  }
  virtual ElementBase& operator=(const ElementBase& other) {
    const Element<T>* const element =  dynamic_cast<const Element<T>*>(&other);
    assert(element != nullptr);
    *this = *element;
    return *this;
  }
  inline virtual int GetDim() const {
    return Traits::kDim;
  }
  virtual void Print() const {
    Traits::Print(GetValue());
  }
  virtual void SetIdentity() {
    Traits::SetIdentity(GetValue());
  }
  virtual void SetRandom(int& s) {
    Traits::SetRandom(GetValue(), s);
  }
  virtual void Boxplus(const VecCRefX& vec, ElementBase* out) const {
    Traits::Boxplus(GetValue(), vec, out->GetValue<T>());
  }
  virtual void Boxminus(const ElementBase& ref, VecRefX vec) const {
    Traits::Boxminus(GetValue(), ref.GetValue<T>(), vec);
  }
  virtual MatX BoxplusJacInp(const VecCRefX& vec) const {
    return Traits::BoxplusJacInp(GetValue(), vec);
  }
  virtual MatX BoxplusJacVec(const VecCRefX& vec) const {
    return Traits::BoxplusJacVec(GetValue(), vec);
  }
  virtual MatX BoxminusJacInp(const ElementBase& ref) const {
    return Traits::BoxminusJacInp(GetValue(), ref.GetValue<T>());
  }
  virtual MatX BoxminusJacRef(const ElementBase& ref) const {
    return Traits::BoxminusJacRef(GetValue(), ref.GetValue<T>());
  }
  T& GetValue() {
    return x_;
  }
  const T& GetValue() const {
    return x_;
  }
 protected:
  T x_;
  const ElementDescription<T>* description_;
};

// ==================== Traits Implementation ==================== //
/*! \brief Scalar Trait.
 *         Element trait for regular scalar.
 */
template<>
class ElementTraits<double> {
 public:
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = 1;
  static void Print(const double& x) {
    std::cout << x << std::endl;
  }
  static const double Identity() {
    return 0;
  }
  static void SetIdentity(double& x) {
    x = 0;
  }
  static void SetRandom(double& x, int& s) {
    static std::default_random_engine generator(s++);
    static std::normal_distribution<double> distribution(0.0, 1.0);
    x = distribution(generator);
  }
  static void Boxplus(const double& in, const VecCRef<kDim>& vec, double& out) {
    out = in + vec(0);
  }
  static void Boxminus(const double& in, const double& ref, VecRef<kDim> vec) {
    vec(0) = in - ref;
  }
  static Mat<kDim> BoxplusJacInp(const double& in, const VecCRef<kDim>& vec) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxplusJacVec(const double& in, const VecCRef<kDim>& vec) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacInp(const double& in, const double& ref) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacRef(const double& in, const double& ref) {
    return -Mat<kDim>::Identity();
  }
};

/*! \brief Vector Trait.
 *         Element trait for vector of scalars
 */
template<int N>
class ElementTraits<Vec<N>> {
 public:
  static constexpr bool kIsVectorSpace = true;
  static constexpr int kDim = N;
  static void Print(const Vec<N>& x) {
    std::cout << x.transpose() << std::endl;
  }
  static const Vec<N> Identity() {
    return Vec<N>::Zero();
  }
  static void SetIdentity(Vec<N>& x) {
    x.setZero();
  }
  static void SetRandom(Vec<N>& x, int& s) {
    static std::default_random_engine generator(s++);
    static std::normal_distribution<double> distribution(0.0, 1.0);
    for (unsigned int i = 0; i < N; i++) {
      x(i) = distribution(generator);
    }
  }
  static void Boxplus(const Vec<N>& in, const VecCRef<kDim>& vec, Vec<N>& out) {
    out = in + vec;
  }
  static void Boxminus(const Vec<N>& in, const Vec<N>& ref, VecRef<kDim> vec) {
    vec = in - ref;
  }
  static Mat<kDim> BoxplusJacInp(const Vec<N>& in, const VecCRef<kDim>& vec) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxplusJacVec(const Vec<N>& in, const VecCRef<kDim>& vec) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacInp(const Vec<N>& in, const Vec<N>& ref) {
    return Mat<kDim>::Identity();
  }
  static Mat<kDim> BoxminusJacRef(const Vec<N>& in, const Vec<N>& ref) {
    return -Mat<kDim>::Identity();
  }
};

/*! \brief Array Trait.
 *         Element trait for array of specific sub-elements
 */
template<typename T, size_t N>
class ElementTraits<std::array<T, N>> {
 public:
  typedef ElementTraits<T> Traits;
  using array = std::array<T, N>;
  static constexpr bool kIsVectorSpace = Traits::kIsVectorSpace;
  static constexpr int kElementDim = Traits::kDim;
  static constexpr int kDim = N * kElementDim;
  static void Print(const array& x) {
    for (const T& i : x) {
      Traits::Print(i);
    }
  }
  static const array Identity() {
    array x;
    SetIdentity(x);
    return x;
  }
  static void SetIdentity(array& x) {
    for (T& i : x) {
      Traits::SetIdentity(i);
    }
  }
  static void SetRandom(array& x, int& s) {
    for (T& i : x) {
      Traits::SetRandom(i, s);
    }
  }
  static void Boxplus(const array& in, const VecCRef<kDim>& vec, array& out) {
    for (int i = 0; i < N; i++) {
      Traits::Boxplus(in[i], vec.template block<kElementDim, 1>(i * kElementDim, 0), out[i]);
    }
  }
  static void Boxminus(const array& in, const array& ref, VecRef<kDim> vec) {
    for (int i = 0; i < N; i++) {
      Traits::Boxminus(in[i], ref[i], vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
  }
  static Mat<kDim> BoxplusJacInp(const array& in, const VecCRef<kDim>& vec) {
    Mat<kDim> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxplusJacInp(in[i], vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
    return J;
  }
  static Mat<kDim> BoxplusJacVec(const array& in, const VecCRef<kDim>& vec) {
    Mat<kDim> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxplusJacVec(in[i], vec.template block<kElementDim, 1>(i * kElementDim, 0));
    }
    return J;
  }
  static Mat<kDim> BoxminusJacInp(const array& in, const array& ref) {
    Mat<kDim> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxminusJacInp(in[i], ref[i]);
    }
    return J;
  }
  static Mat<kDim> BoxminusJacRef(const array& in, const array& ref) {
    Mat<kDim> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<kElementDim, kElementDim>(i * kElementDim, i * kElementDim) =
          Traits::BoxminusJacRef(in[i], ref[i]);
    }
    return J;
  }
};

/*! \brief Unit Quaternion Trait.
 *         Element trait for unit quaternion. Employed to represent
 *         orientations.
 */
template<>
class ElementTraits<Quat> {
 public:
  static constexpr bool kIsVectorSpace = false;
  static constexpr int kDim = 3;
  static void Print(const Quat& x) {
    std::cout << x << std::endl;
  }
  static const Quat Identity() {
    return Quat();
  }
  static void SetIdentity(Quat& x) {
    x.setIdentity();
  }
  static void SetRandom(Quat& x, int& s) {
    static std::default_random_engine generator(s++);
    static std::normal_distribution<double> distribution(0.0, 1.0);
    x.toImplementation().w() = distribution(generator);
    x.toImplementation().x() = distribution(generator);
    x.toImplementation().y() = distribution(generator);
    x.toImplementation().z() = distribution(generator);
    x.fix();
  }
  static void Boxplus(const Quat& in, const VecCRef<kDim>& vec, Quat& out) {
    out = in.boxPlus(vec);
  }
  static void Boxminus(const Quat& in, const Quat& ref, VecRef<kDim> vec) {
    vec = in.boxMinus(ref);
  }
  static Mat<kDim> BoxplusJacInp(const Quat& in, const VecCRef<kDim>& vec) {
    RotMat m = m.exponentialMap(vec);
    return m.matrix();
  }
  static Mat<kDim> BoxplusJacVec(const Quat& in, const VecCRef<kDim>& vec) {
    return GammaMat(vec);
  }
  static Mat<kDim> BoxminusJacInp(const Quat& in, const Quat& ref) {
    return GammaMat(in.boxMinus(ref)).inverse();
  }
  static Mat<kDim> BoxminusJacRef(const Quat& in, const Quat& ref) {
    return -GammaMat(in.boxMinus(ref)).inverse()* RotMat(in * ref.inverted()).matrix();
  }
};

}

#endif /* GIF_ELEMENT_HPP_ */
