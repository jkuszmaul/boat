#include "dynamics.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

DEFINE_double(kp, 4.0, "kP");
DEFINE_double(kd, 2.0, "kD");
DEFINE_double(kpm, 10.0, "kPM");
DEFINE_double(kdm, 1.0, "kDM");

namespace aero {

using Eigen::Quaterniond;

void MotorVels(double thrust, const Array3d &moments, const Param &p,
               Input *motor_vels) {
  if (motor_vels == nullptr) {
    return;
  }
  double gamma = p.kM / p.kF;
  Eigen::Vector4d outs;
  outs << thrust, moments.x(), moments.y(), moments.z();
  Eigen::Matrix4d F2Inputs;
  F2Inputs << 1, 1, 1, 1,
              0, p.L, 0, -p.L,
              -p.L, 0, p.L, 0,
              gamma, -gamma, gamma, -gamma;
  Eigen::Vector4d motor_forces = F2Inputs.inverse() * outs;
  *motor_vels = motor_forces / p.kF;
  *motor_vels = motor_vels->sqrt();
}

void ForcesMoments(Input motor_vels, const Param &p, Array3d *forces,
                   Array3d *moments) {
  if (forces == nullptr || moments == nullptr) {
    return;
  }
  // Impose saturation limits:
  motor_vels = motor_vels.max(p.minvin);
  motor_vels = motor_vels.min(p.maxvin);

  // Compute thrust
  Input motor_forces = motor_vels.square() * p.kF;
  *forces << 0,
             0,
             -motor_forces.sum();

  // Compute moments
  double gamma = p.kM / p.kF;
  Eigen::Matrix<double, 3, 4> F2M;
  F2M << 0.0, p.L, 0.0, -p.L,
         -p.L, 0.0, p.L, 0.0,
         gamma, -gamma, gamma, -gamma;
  *moments = F2M * motor_forces.matrix();
}

void Dynamics(const State &x, const Array3d &forces, const Array3d &moments, const Param &p, State *xdot) {
  if (xdot == nullptr) {
    return;
  }

  Eigen::Matrix<double, 4, 1> qmat = x.block(6, 0, 4, 1);
  Quaterniond q(qmat);
  q.normalize();
  Vector3d Omega = x.block(10, 0, 3, 1);
  Array3d Fi =
      (q * Quaterniond(0, forces.x(), forces.y(), forces.z()) * q.inverse()).vec();
  // Derivative of position is velocity
  xdot->block(0, 0, 3, 1) = x.block(3, 0, 3, 1);
  // Accelerations:
  xdot->block(3, 0, 3, 1) = Fi / p.m + Array3d(0, 0, p.g);
  // Derivative of quaternion
  xdot->block(6, 0, 4, 1) =
      0.5 * (q * Quaterniond(0, Omega.x(), Omega.y(), Omega.z())).coeffs();
  // Angular Accelerations
  xdot->block(10, 0, 3, 1) =
      p.J.inverse() * (moments.matrix() - Omega.cross(p.J * Omega));
}

}  // namespace aero
