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
  static constexpr bool is_vectorspace_ = true;
  static constexpr int d_ = 0;
  static void print(const T& x) {
  }
  static const T identity() {
    T x;
    return x;
  }
  static void setIdentity(T& x) {
  }
  static void setRandom(T& x, int& s) {
  }
  static void boxplus(const T& ref, const VecRC<d_>& vec, T& out) {
    out = ref;   // Must be computable in-place
  }
  static void boxminus(const T& in, const T& ref, VecR<d_> vec) {
  }
  static Mat<d_, d_> boxplusJacInp(const T& in, const VecRC<d_>& vec) {
    return Mat<d_, d_>::Identity();
  }
  static Mat<d_, d_> boxplusJacVec(const T& in, const VecRC<d_>& vec) {
    return Mat<d_, d_>::Identity();
  }
  static Mat<d_, d_> boxminusJacInp(const T& in, const T& ref) {
    return Mat<d_, d_>::Identity();
  }
  static Mat<d_, d_> boxminusJacRef(const T& in, const T& ref) {
    return Mat<d_, d_>::Identity();
  }
};

/*! \brief Element Base.
 *         Will be used to store a vector of elements in ElementVector
 */
class ElementBase {
 public:
  ElementBase() {
  }
  virtual ~ElementBase() {
  }
  virtual ElementBase& operator=(const ElementBase& other) = 0;
  virtual int getDim() const = 0;
  virtual void print() const = 0;
  virtual void setIdentity() = 0;
  virtual void setRandom(int& s) = 0;
  virtual void boxplus(const VecRC<>& vec,const SP<ElementBase>& out) const = 0;
  virtual void boxminus(const SP<const ElementBase>& ref, VecR<> vec) const = 0;
  virtual Mat<> boxplusJacInp(const VecRC<>& vec) const = 0;
  virtual Mat<> boxplusJacVec(const VecRC<>& vec) const = 0;
  virtual Mat<> boxminusJacInp(const SP<ElementBase>& ref) const = 0;
  virtual Mat<> boxminusJacRef(const SP<ElementBase>& ref) const = 0;
  template<typename T>
  T& get() {
    return dynamic_cast<Element<T>*>(this)->get();
  }
  template<typename T>
  const T& get() const {
    return dynamic_cast<const Element<T>*>(this)->get();
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
  }
  virtual ~Element() {
  }
  virtual Element<T>& operator=(const Element<T>& other) {
    get() = other.get();
    return *this;
  }
  virtual ElementBase& operator=(const ElementBase& other) {
    *this = dynamic_cast<const Element<T>&>(other);
    return *this;
  }
  inline virtual int getDim() const {
    return Traits::d_;
  }
  virtual void print() const {
    Traits::print(get());
  }
  virtual void setIdentity() {
    Traits::setIdentity(get());
  }
  virtual void setRandom(int& s) {
    Traits::setRandom(get(), s);
  }
  virtual void boxplus(const VecRC<>& vec, const SP<ElementBase>& out) const {
    Traits::boxplus(get(), vec, out->get<T>());
  }
  virtual void boxminus(const SP<const ElementBase>& ref, VecR<> vec) const {
    Traits::boxminus(get(), ref->get<T>(), vec);
  }
  virtual Mat<> boxplusJacInp(const VecRC<>& vec) const {
    return Traits::boxplusJacInp(get(), vec);
  }
  virtual Mat<> boxplusJacVec(const VecRC<>& vec) const {
    return Traits::boxplusJacVec(get(), vec);
  }
  virtual Mat<> boxminusJacInp(const SP<ElementBase>& ref) const {
    return Traits::boxminusJacInp(get(), ref->get<T>());
  }
  virtual Mat<> boxminusJacRef(const SP<ElementBase>& ref) const {
    return Traits::boxminusJacRef(get(), ref->get<T>());
  }
  T& get() {
    return x_;
  }
  const T& get() const {
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
  static constexpr bool is_vectorspace_ = true;
  static constexpr int d_ = 1;
  static void print(const double& x) {
    std::cout << x << std::endl;
  }
  static const double identity() {
    return 0;
  }
  static void setIdentity(double& x) {
    x = 0;
  }
  static void setRandom(double& x, int& s) {
    std::default_random_engine generator(s++);
    std::normal_distribution<double> distribution(0.0, 1.0);
    x = distribution(generator);
  }
  static void boxplus(const double& in, const VecRC<d_>& vec, double& out) {
    out = in + vec(0);
  }
  static void boxminus(const double& in, const double& ref, VecR<d_> vec) {
    vec(0) = in - ref;
  }
  static Mat<d_> boxplusJacInp(const double& in, const VecRC<d_>& vec) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxplusJacVec(const double& in, const VecRC<d_>& vec) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxminusJacInp(const double& in, const double& ref) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxminusJacRef(const double& in, const double& ref) {
    return -Mat<d_>::Identity();
  }
};

/*! \brief Vector Trait.
 *         Element trait for vector of scalars
 */
template<int N>
class ElementTraits<Vec<N>> {
 public:
  static constexpr bool is_vectorspace_ = true;
  static constexpr int d_ = N;
  static void print(const Vec<N>& x) {
    std::cout << x.transpose() << std::endl;
  }
  static const Vec<N> identity() {
    return Vec<N>::Zero();
  }
  static void setIdentity(Vec<N>& x) {
    x.setZero();
  }
  static void setRandom(Vec<N>& x, int& s) {
    std::default_random_engine generator(s++);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (unsigned int i = 0; i < N; i++) {
      x(i) = distribution(generator);
    }
  }
  static void boxplus(const Vec<N>& in, const VecRC<d_>& vec, Vec<N>& out) {
    out = in + vec;
  }
  static void boxminus(const Vec<N>& in, const Vec<N>& ref, VecR<d_> vec) {
    vec = in - ref;
  }
  static Mat<d_> boxplusJacInp(const Vec<N>& in, const VecRC<d_>& vec) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxplusJacVec(const Vec<N>& in, const VecRC<d_>& vec) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxminusJacInp(const Vec<N>& in, const Vec<N>& ref) {
    return Mat<d_>::Identity();
  }
  static Mat<d_> boxminusJacRef(const Vec<N>& in, const Vec<N>& ref) {
    return -Mat<d_>::Identity();
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
  static constexpr bool is_vectorspace_ = Traits::is_vectorspace_;
  static constexpr int ed_ = Traits::d_;
  static constexpr int d_ = N * ed_;
  static void print(const array& x) {
    for (const T& i : x) {
      Traits::print(i);
    }
  }
  static const array identity() {
    array x;
    setIdentity(x);
    return x;
  }
  static void setIdentity(array& x) {
    for (T& i : x) {
      Traits::setIdentity(i);
    }
  }
  static void setRandom(array& x, int& s) {
    for (T& i : x) {
      Traits::setRandom(i, s);
    }
  }
  static void boxplus(const array& in, const VecRC<d_>& vec, array& out) {
    for (int i = 0; i < N; i++) {
      Traits::boxplus(in[i], vec.template block<ed_, 1>(i * ed_, 0), out[i]);
    }
  }
  static void boxminus(const array& in, const array& ref, VecR<d_> vec) {
    for (int i = 0; i < N; i++) {
      Traits::boxminus(in[i], ref[i], vec.template block<ed_, 1>(i * ed_, 0));
    }
  }
  static Mat<d_> boxplusJacInp(const array& in, const VecRC<d_>& vec) {
    Mat<d_> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<ed_, ed_>(i * ed_, i * ed_) =
          Traits::boxplusJacInp(in[i], vec.template block<ed_, 1>(i * ed_, 0));
    }
    return J;
  }
  static Mat<d_> boxplusJacVec(const array& in, const VecRC<d_>& vec) {
    Mat<d_> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<ed_, ed_>(i * ed_, i * ed_) =
          Traits::boxplusJacVec(in[i], vec.template block<ed_, 1>(i * ed_, 0));
    }
    return J;
  }
  static Mat<d_> boxminusJacInp(const array& in, const array& ref) {
    Mat<d_> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<ed_, ed_>(i * ed_, i * ed_) =
          Traits::boxminusJacInp(in[i], ref[i]);
    }
    return J;
  }
  static Mat<d_> boxminusJacRef(const array& in, const array& ref) {
    Mat<d_> J;
    J.setZero();
    for (int i = 0; i < N; i++) {
      J.template block<ed_, ed_>(i * ed_, i * ed_) =
          Traits::boxminusJacRef(in[i], ref[i]);
    }
    return J;
  }
};

/*! \brief Unit Quaternion Trait.
 *         Element trait for unit quaternion. Employed to represent
 *         orientations.
 */
template<>
class ElementTraits<QPD> {
 public:
  static constexpr bool is_vectorspace_ = false;
  static constexpr int d_ = 3;
  static void print(const QPD& x) {
    std::cout << x << std::endl;
  }
  static const QPD identity() {
    return QPD();
  }
  static void setIdentity(QPD& x) {
    x.setIdentity();
  }
  static void setRandom(QPD& x, int& s) {
    std::default_random_engine generator(s++);
    std::normal_distribution<double> distribution(0.0, 1.0);
    x.toImplementation().w() = distribution(generator);
    x.toImplementation().x() = distribution(generator);
    x.toImplementation().y() = distribution(generator);
    x.toImplementation().z() = distribution(generator);
    x.fix();
  }
  static void boxplus(const QPD& in, const VecRC<d_>& vec, QPD& out) {
    out = in.boxPlus(vec);
  }
  static void boxminus(const QPD& in, const QPD& ref, VecR<d_> vec) {
    vec = in.boxMinus(ref);
  }
  static Mat<d_> boxplusJacInp(const QPD& in, const VecRC<d_>& vec) {
    MPD m = m.exponentialMap(vec);
    return m.matrix();
  }
  static Mat<d_> boxplusJacVec(const QPD& in, const VecRC<d_>& vec) {
    return Lmat(vec);
  }
  static Mat<d_> boxminusJacInp(const QPD& in, const QPD& ref) {
    return Lmat(in.boxMinus(ref)).inverse();
  }
  static Mat<d_> boxminusJacRef(const QPD& in, const QPD& ref) {
    return -Lmat(in.boxMinus(ref)).inverse()* MPD(in * ref.inverted()).matrix();
  }
};

}

#endif /* GIF_ELEMENT_HPP_ */
