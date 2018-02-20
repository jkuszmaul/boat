#pragma once

#include <Eigen/Core>

#include "params.h"

namespace aero {

typedef Eigen::Array<double, 13, 1> State;
typedef Eigen::Array<double, 4, 1> Input;
typedef Eigen::Array<double, 3, 1> Array3d;
typedef Eigen::Vector3d Vector3d;

// Does the reverse of ForcesMoments, computing the needed motor velocities to
// achieve a given thrust and set of moments.
void MotorVels(double thrust, const Array3d &moments, const Param &p,
               Input *motor_vels);

// Compute the body-frame forces and moments of the quadrotor given the input
// motor velocities.
void ForcesMoments(Input motor_vels, const Param &p, Array3d *forces,
                   Array3d *moments);

// A function which, given the body-frame forces and moments on
// a given system, excluding gravity (we will add in gravity ourselves),
// computes the derivative of the system state.
// The system state shall be a column vector with the following
// elements, in order:
// 1) The inertial frame position (x, y, z)
// 2) The inertial frame velocity (vx, vy, vz)
// 3) The attitude quaternion (x, y, z, w), with w = scalar component (because Eigen is stupid)
// 4) The body-frame rotational velocities (omega_x, omega_y, omega_z)
void Dynamics(const State &x, const Array3d &forces, const Array3d &moments,
              const Param &p, State *xdot);

}  // namespace aero
