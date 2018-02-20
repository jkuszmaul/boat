#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "boost/numeric/odeint.hpp"
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include "boost/numeric/odeint/external/eigen/eigen_resize.hpp"

namespace boost {
namespace numeric {
namespace odeint {
template <typename B, int S1, int S2, int O, int M1, int M2>
struct vector_space_norm_inf<Eigen::Array<B, S1, S2, O, M1, M2>> {
  typedef B result_type;
  result_type operator()(const Eigen::Array<B, S1, S2, O, M1, M2> &m) const {
    return Eigen::abs(m).maxCoeff();
  }
};

template <class Derived>
struct algebra_dispatcher_sfinae<
    Derived, typename boost::enable_if<typename boost::is_base_of<
                 Eigen::MatrixBase<Derived>, Derived>::type>::type> {
  typedef vector_space_algebra algebra_type;
};

template <class Derived>
struct algebra_dispatcher_sfinae<
    Derived, typename boost::enable_if<typename boost::is_base_of<
                 Eigen::ArrayBase<Derived>, Derived>::type>::type> {
  typedef vector_space_algebra algebra_type;
};
}
}
}
