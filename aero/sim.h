#pragma once

#include "control.h"
#include "dynamics.h"
#include "eigen_int.hpp"
#include "gflags.h"
#include <functional>
#include <map>

//DEFINE_double(thrust, 0.0, "Thrust to give to dynamics");
//DEFINE_double(mx, 0.0, "X moment");
//DEFINE_double(my, 0.0, "Y moment");
//DEFINE_double(mz, 0.0, "z moment");
//DEFINE_double(m1, 0.0, "Motor 1");
//DEFINE_double(m2, 0.0, "Motor 2");
//DEFINE_double(m3, 0.0, "Motor 3");
//DEFINE_double(m4, 0.0, "Motor 4");

namespace aero {
using namespace boost::numeric::odeint;

using Eigen::Quaterniond;
class Simulate {
 public:
  Simulate() :
    //traj_({0, .1},
    //      {0},
    //      {-16.764, 0.1}) {
    traj_({0, 14.469, -1.0127e-1, -1.1841e-4, -7.6013e-4, -3.4094e-3, 2.9767e-4},
          {9.144, 3.8769, -4.1309, 1.1727, -1.6510e-1, 1.2217e-2, -3.828e-4},
          {-16.764, 7.8507e-1, 3.1597, -1.0059, 1.2597e-1, -6.9039e-3, 1.2789e-4}) {
  }

  void Run(double tf) {
    double niter = 1500;
    double dt = tf / niter;
    State x0 = CalcX0();
    //std::cout << "x0: " << x0.transpose() << "\n";

    //dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<State>>>
    //    integrator;
    runge_kutta_dopri5<State> dopri5;;
    //integrator.initialize(x0, 0, tf / niter);
    auto f = std::bind(&Simulate::Step, this, std::placeholders::_1,
                       std::placeholders::_2, std::placeholders::_3);
    State xint = x0;
    double t = 0;
    std::cout << "#t,px,py,pz,vx,vy,vz,qx,qy,qz,qw,wx,wy,wz,accx,accy,accz,trajx,trajy,trajz,trajvx,trajvy,trajvz,Fbz,mx,my,mz,qdx,qdy,qdz,qdw";
    State lastx;
    for (double ii = 0; ii <= niter; ii += 1) {
    //while (t < tf) {
      //std::pair<double, double> tstep = integrator.do_step(f);
      //t = tstep.second;
      if (true) {
        lastx = xint;
        dopri5.do_step(f, xint, t, dt);
      } else {
        State xdot;
        lastx = xint;
        Step(xint, xdot, t);
        xint += xdot * dt;
      }
      State x = lastx;
      State xdot = xdots_.lower_bound(t)->second;
      //std::cout << "t0: " << tstep.first << " tf: " << tstep.second << std::endl;
      //integrator.calc_state(t, x);
      for (int jj = 0; jj < 13; ++jj) {
        if (std::abs(x(jj, 0)) < 1e-5) {
          x(jj, 0) = 0;
        }
        if (std::abs(xdot(jj, 0)) < 1e-8) {
          xdot(jj, 0) = 0;
        }
      }
      x.block(6, 0, 4, 1) /=
          std::sqrt(x[6] * x[6] + x[7] * x[7] + x[8] * x[8] + x[9] * x[9]);
      //std::cout << t << " xdot: " << xdot.transpose() << std::endl;
      Array3d pos, vel, accel;
      traj_.PolyPos(t, 0, &pos);
      traj_.PolyPos(t, 1, &vel);
      traj_.PolyPos(t, 2, &accel);
      auto it = goal_forces.lower_bound(t);
      auto elem = it->second;
      auto fm = forces_moments.lower_bound(t)->second;
      std::cout << t << " " << lastx.transpose() << " "
                << xdot.block(3, 0, 3, 1).transpose() << " " << pos.transpose()
                << " " << vel.transpose() << " " << fm.first[2] << " "
                << fm.second.transpose() << " "
                << qds_.lower_bound(t)->second.coeffs().transpose() << "\n";
      //std::cout << t << " traj: " << pos.transpose()
      //          << " traj vel: " << vel.transpose()
      //          << " goal time: " << it->first
      //          << " thrust: " << elem.first
      //          << " moment: " << elem.second.transpose()
      //          << " true force: " << fm.first.transpose()
      //          << " true moment: " << fm.second.transpose()
      //          << std::endl;
      //std::cout << t << ": " << x.transpose() << std::endl;
      t += dt;
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
        Eigen::Quaterniond(std::cos(psi / 2.0), 0, 0, std::sin(-psi / 2.0)).coeffs();
    x0.block(6, 0, 4, 1) << -.55, .09, .1, .82;
    x0.block(10, 0, 3, 1) = Array3d(0, 0, -psidot);
    //x0 *= 0;
    //x0(9, 0) = 1.0;
    return x0;
  }
  void KalmanStep(const Eigen::Array<double, NSTATE, NSTATE + 1> &PX,
                  Eigen::Array<double, NSTATE, NSTATE + 1> &PXdot, double t) {
    State x = PX.col(NSTATE);
    Step(x, PXdot.col(NSTATE), t);
  }
  void Step(const State &x, State &xdot, double t) {
    Array3d pos_d, vel_d, accel_d, jerk_d;
    double psi_d, psidot_d;
    traj_.StateAtTime(t, &pos_d, &vel_d, &accel_d, &jerk_d, &psi_d, &psidot_d);

    Array3d moment_des;
    double thrust_des;
    Quaterniond qd_out;
    GoalForces(x, pos_d, vel_d, accel_d, jerk_d, psi_d, psidot_d, p_, &thrust_des, &moment_des, &qd_out);
    goal_forces[t] = {thrust_des, moment_des};
    qds_[t] = qd_out;

    //thrust_des = FLAGS_thrust;
    //moment_des << FLAGS_mx, FLAGS_my, FLAGS_mz;

    Input motor_velocities;
    MotorVels(thrust_des, moment_des, p_, &motor_velocities);

    //motor_velocities << FLAGS_m1, FLAGS_m2, FLAGS_m3, FLAGS_m4;

    Array3d F, M;
    ForcesMoments(motor_velocities, p_, &F, &M);
    forces_moments[t] = {F, M};

    //F << 0, 0, FLAGS_thrust;
    //M << FLAGS_mx, FLAGS_my, FLAGS_mz;
    Dynamics(x, F, M, p_, &xdot);
    xdots_[t] = State(xdot);
  }
  Trajectory traj_;
  Param p_;
  Sensors sensors_;

  std::map<double, std::pair<double, Array3d>> goal_forces;
  std::map<double, std::pair<Array3d, Array3d>> forces_moments;
  std::map<double, State> xdots_;
  std::map<double, Quaterniond> qds_;
};  // class Simulate
}  // namespace aero
