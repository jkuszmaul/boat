#pragma once

#include "control.h"
#include "dynamics.h"
#include "eigen_int.hpp"
#include "gflags.h"
#include <functional>

DEFINE_double(thrust, 0.0, "Thrust to give to dynamics");
DEFINE_double(mx, 0.0, "X moment");
DEFINE_double(my, 0.0, "Y moment");
DEFINE_double(mz, 0.0, "z moment");

namespace aero {
using namespace boost::numeric::odeint;

class Simulate {
 public:
  Simulate() :
    traj_({0, 14.469, -1.0127e-1, -1.1841e-4, -7.6013e-4, -3.4094e-3, 2.9767e-4},
          {9.144, 3.8769, -4.1309, 1.1727, -1.6510e-1, 1.2217e-2, -3.828e-4},
          {-16.764, 7.8507e-1, 3.1597, -1.0059, 1.2597e-1, -6.9039e-3, 1.2789e-4}) {
  }

  void Run(double dt) {
    State x0 = CalcX0();

    dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<State>>>
        integrator;
    p_.g = 0;
    integrator.initialize(x0, 0, dt);
    std::pair<double, double> tstep = integrator.do_step(
        std::bind(&Simulate::Step, this, std::placeholders::_1,
                  std::placeholders::_2, std::placeholders::_3));
    std::cout << "t0: " << tstep.first << " tf: " << tstep.second << std::endl;
    State x;
    for (double ii = 0; ii <= 1; ii += 0.001) {
      double t = ii * dt;
      integrator.calc_state(t, x);
      for (int jj = 0; jj < 13; ++jj) {
        if (std::abs(x(jj, 0)) < 1e-5) {
          x(jj, 0) = 0;
        }
      }
      x.block(6, 0, 4, 1) /=
          std::sqrt(x[6] * x[6] + x[7] * x[7] + x[8] * x[8] + x[9] * x[9]);
      std::cout << t << ": " << x.transpose() << std::endl;
//      Array3d pos;
//      traj_.PolyPos(t, 0, &pos);
//      std::cout << t << " traj: " << pos.transpose() << std::endl;
    }
  }

 private:
  State CalcX0() {
    Array3d pos, vel;
    double psi, psidot;
    traj_.StateAtTime(0, &pos, &vel, nullptr, nullptr, &psi, &psidot);
    State x0;
    x0.block(0, 0, 3, 1) = pos;
    x0.block(3, 0, 3, 1) = vel;
    // Note that the Z-axis points down, and so psi = -rot(z)
    x0.block(6, 0, 4, 1) =
        Eigen::Quaterniond(std::cos(psi / 2.0), 0, 0, std::sin(psi / 2.0)).coeffs();
    x0.block(10, 0, 3, 1) = Array3d(0, 0, -psidot);
    x0 *= 0;
    x0(9, 0) = 1.0;
    return x0;
  }
  void Step(const State &x, State &xdot, double t) {
    Array3d pos_d, vel_d, accel_d, jerk_d;
    double psi_d, psidot_d;
    traj_.StateAtTime(t, &pos_d, &vel_d, &accel_d, &jerk_d, &psi_d, &psidot_d);

    Array3d moment_des;
    double thrust_des;
    GoalForces(x, pos_d, vel_d, accel_d, jerk_d, psi_d, psidot_d, p_, &thrust_des, &moment_des);

    Input motor_velocities;
    MotorVels(thrust_des, moment_des, p_, &motor_velocities);

    Array3d F, M;
    ForcesMoments(motor_velocities, p_, &F, &M);

    F << 0, 0, FLAGS_thrust;
    M << FLAGS_mx, FLAGS_my, FLAGS_mz;
    Dynamics(x, F, M, p_, &xdot);
  }
  Trajectory traj_;
  Param p_;
};  // class Simulate
}  // namespace aero
