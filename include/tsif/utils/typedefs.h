#ifndef TSIF_TYPEDEFS_HPP_
#define TSIF_TYPEDEFS_HPP_

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace tsif{

template<int N = -1>
using Vec = Eigen::Matrix<double,N,1>;
using Vec2 = Vec<2>;
using Vec3 = Vec<3>;
using VecX = Vec<>;
template<int N = -1>
using VecRef = Eigen::Ref<Vec<N>>;
using VecRef2 = VecRef<2>;
using VecRef3 = VecRef<3>;
using VecRefX = VecRef<>;
template<int N = -1>
using VecCRef = Eigen::Ref<const Vec<N>>;
using VecCRef2 = VecCRef<2>;
using VecCRef3 = VecCRef<3>;
using VecCRefX = VecCRef<>;

template<int N = -1, int M = N>
using Mat = Eigen::Matrix<double,N,M>;
using Mat2 = Mat<2>;
using Mat3 = Mat<3>;
using MatX = Mat<>;
template<int N = -1, int M = N>
using MatRef = typename std::conditional<(N==1 & M>1),
                                         Eigen::Ref<Mat<N,M>,0,Eigen::InnerStride<>>,
                                         Eigen::Ref<Mat<N,M>>>::type;
using MatRef2 = MatRef<2>;
using MatRef3 = MatRef<3>;
using MatRefX = MatRef<>;
template<int N = -1, int M = N>
using MatCRef = Eigen::Ref<const Mat<N,M>>;
using MatCRef2 = MatCRef<2>;
using MatCRef3 = MatCRef<3>;
using MatCRefX = MatCRef<>;

typedef Eigen::Quaterniond Quat;

} // namespace tsif

#endif /* TSIF_TYPEDEFS_HPP_ */
