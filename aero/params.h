#pragma once

#include <Eigen/Core>
#include "gflags.h"

DECLARE_double(kp);
DECLARE_double(kd);
DECLARE_double(kpm);
DECLARE_double(kdm);

namespace aero {
struct Param {
  Param()
      : m(0.5), L(0.175), kM(1.5e-9 * 3600.0 / M_PI / M_PI),
        kF(6.11e-8 * 3600.0 / M_PI / M_PI) {
    J = Eigen::Matrix3d::Zero();
    J.diagonal() << 2.32e-3, 2.32e-3, 4e-3;

    kP = Eigen::Matrix3d::Zero();
    kP.diagonal() << 1, 1, 1;
    kP *= FLAGS_kp;

    kD = Eigen::Matrix3d::Zero();
    kD.diagonal() << 1, 1, 1;
    kD *= FLAGS_kd;

    kP_M = Eigen::Matrix3d::Zero();
    kP_M.diagonal() << 1, 1, 1;
    kP_M *= FLAGS_kpm;

    kD_M = Eigen::Matrix3d::Zero();
    kD_M.diagonal() << 1, 1, 1;
    kD_M *= FLAGS_kdm;
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
