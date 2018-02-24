#pragma once
#include <Eigen/Core>
#include <functional>

namespace math {

template <typename T, int N, int M>
Eigen::Matrix<T, M, N>
    jacobian(std::function<Eigen::Matrix<T, M, 1>(Eigen::Matrix<T, N, 1>)> f,
             Eigen::Matrix<T, N, 1> x0) {
  Eigen::Matrix<T, M, 1> y1, y2;
  Eigen::Matrix<T, N, 1> xdiff;
  Eigen::Matrix<T, M, N> J;
  T h = 1e-5;
  for (int ii = 0; ii < N; ++ii) {
    xdiff = xdiff.Zero();
    xdiff(ii, 0) = h;
    y1 = f(x0 + xdiff);
    y2 = f(x0 - xdiff);
    J.col(ii) = (y1 - y2) / (2.0 * h);
  }
  return J;
}

}  // namespace math
