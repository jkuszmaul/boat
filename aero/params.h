#pragma once

#include <Eigen/Core>

namespace aero {
struct Param {
  Param() : m(0.5), L(0.175), kM(1.5e-9), kF(6.11e-8) {
    J.diagonal() << 2.32e-3, 2.32e-3, 4e-3;
    kP.diagonal() << 1, 1, 1;
    kD.diagonal() << 1, 1, 1;
    kP_M.diagonal() << 1, 1, 1;
    kD_M.diagonal() << 1, 1, 1;
  }
  // Acceleration due to gravity
  double g = 9.8;
  // mass
  double m;
  // Moments of Inertia
  Eigen::Matrix3d J;
  // Radius of quadrotor
  double L;
  // Motor constants
  double kM, kF;

  // Input saturation limits
  double minvin = 0, maxvin = 9000 * M_PI / 60;

  // Controller gains
  Eigen::Matrix3d kP, kD;
  Eigen::Matrix3d kP_M, kD_M; // For moments
};  // Param
}  // namespace aero
