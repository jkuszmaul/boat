#include "control.h"
#include "params.h"

#include <cmath>
#include <iostream>

namespace aero {
using Eigen::Quaterniond;
namespace {
  template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
  }
}
void GoalForces(const State &x, const Array3d &r_d, const Array3d &rdot_d,
                const Array3d &rddot_d, const Array3d &rdddot_d,
                const double psi_d, const double psidot_d, const Param &p,
                double *thrust_des, Array3d *moment_des, Quaterniond *qd_out) {
  if (thrust_des == nullptr || moment_des == nullptr) {
    return;
  }
  const Vector3d r = x.block(0, 0, 3, 1); // Inertial frame Position
  const Vector3d rdot = x.block(3, 0, 3, 1); // Inertial frame velocity
  Eigen::Matrix<double, 4, 1> qmat = x.block(6, 0, 4, 1); // Quaternion, Inertial->Body
  Quaterniond q(qmat);
  Vector3d Omega = x.block(10, 0, 3, 1); // Body-frame angular velocities

  Vector3d e = r - r_d.matrix(); // Positional error to reference trajectory
  Vector3d edot = rdot - rdot_d.matrix(); // Velocity error to reference trajectory

  Vector3d rfbddot = -p.kP * e - p.kD * edot; // Feedback term, in acceleration
  Vector3d rddot = rddot_d.matrix() + rfbddot; // Total desired acceleration

  Vector3d eddot = rddot - rddot_d.matrix(); // Second derivative of error
  Vector3d rfbdddot = -p.kP * edot - p.kD * eddot; // Derivative of feedback term (jerk)

  Vector3d Fi = (rddot - Vector3d(0, 0, p.g)) * p.m; // Total external (non-gravitational) forces, inertial frame
  *thrust_des = Fi.norm(); // Total needed thrust

  Vector3d Fbbar(0, 0, -1); // Normalized body-frame thrust vector
  Vector3d Fibar = Fi.normalized(); // Normalized inertial-frame thrust vector
  // Compute orientation needed to achieve desired thrust:
  double Fbi = Fbbar.dot(Fibar);
  double qdscale = 1.0 / std::sqrt(2.0 * (1.0 + Fbi));
  Quaterniond qdtilde;
  qdtilde.w() = (1.0 + Fbi) * qdscale;
  qdtilde.vec() = Fbbar.cross(Fibar) * qdscale;
  // Rotate orientation to account for desired heading angle
  Quaterniond qd =
      qdtilde * Quaterniond(std::cos(psi_d / 2.0), 0, 0, std::sin(-psi_d / 2.0));
  //std::cout << "Fi: " << Fi.transpose() << " qd: " << qd.coeffs().transpose() << "\n";
  // Derivative of inertial frame force, for calculating desired angular rotations
  Vector3d Fidot = p.m * (rdddot_d.matrix() + rfbdddot);
  // Compute desired angular rotation
  double fnorm = Fi.norm();
  Vector3d Fibardot = Fidot / fnorm - Fi * Fi.dot(Fidot) / (fnorm * fnorm * fnorm);
  Vector3d Fcross = Fibar.cross(Fibardot);
  Vector3d Omega_des(Fcross.x(), Fcross.y(), psidot_d);
  //Omega_des *= 0;
  // Error quaternion:
  Quaterniond qe = q.conjugate() * qd;
  if (qd_out != nullptr) {
    *qd_out = qd;
  }
  // Feedback term on moments:
  *moment_des = (sgn(qe.w()) * p.kP_M * qe.vec() - p.kD_M * (Omega - Omega_des)).array();
}


}  // namespace aero
