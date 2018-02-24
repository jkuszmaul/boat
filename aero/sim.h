#pragma once

#include "control.h"
#include "dynamics.h"
#include "eigen_int.hpp"
#include "gflags.h"
#include <functional>
#include "math/jacobian.h"
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
    //double niter = 4;
    double dt = tf / niter;
    int ngps = 0.3 / dt, nimu = 0.02 / dt;
    State x0 = CalcX0();
    // contains P and xhat
    Eigen::Array<double, NSTATE, NSTATE+1> PX, PXdot;
    PX = PX.Zero();
    for (int ii = 0; ii < NSTATE; ++ii) {
      PX(ii, ii) = 0.1;
    }
    PX.col(NSTATE) = x0;
    PX.block(0, NSTATE, 3, 1) += 0.1;
    PX.block(3, NSTATE, 3, 1) -= 0.01;
    PX.block(6, NSTATE, 4, 1) -= 0;
    PX.block(10, NSTATE, 3, 1) *= 1.05;
    //std::cout << "x0: " << x0.transpose() << "\n";

    //dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<State>>>
    //    integrator;
    runge_kutta_dopri5<State> dopri5;;
    //integrator.initialize(x0, 0, tf / niter);
    auto f = std::bind(&Simulate::Step, this, std::placeholders::_1,
                       std::placeholders::_1, std::placeholders::_2,
                       std::placeholders::_3);
    State xint = x0;
    double t = 0;
    std::cout << "#t,px,py,pz,vx,vy,vz,qx,qy,qz,qw,wx,wy,wz,accx,accy,accz,"
                 "trajx,trajy,trajz,trajvx,trajvy,trajvz,Fbz,mx,my,mz,qdx,qdy,"
                 "qdz,qdw,pxhat,pyhat,pzhat,vxhat,vyhat,vzhat,qxhat,qyhat,"
                 "qzhat,qwhat,wxhat,wyhat,wzhat\n";
    State lastx;
    for (double ii = 0; ii <= niter; ii += 1) {
    //while (t < tf) {
      //std::pair<double, double> tstep = integrator.do_step(f);
      //t = tstep.second;
      if (false) {
        lastx = xint;
        dopri5.do_step(f, xint, t, dt);
      } else {
        State xdot;
        lastx = xint;
        Step(xint, PX.col(NSTATE), xdot, t);
        //Step(xint, xint, xdot, t);
        //std::cout << xdot.transpose() << "\n";
        xint += xdot * dt;

        KalmanPredict(PX, PXdot, t);
        //std::cout << PXdot.col(NSTATE).transpose() << "\n";
        PX += PXdot * dt;
        //std::cout << xint.transpose() << "\n";
        //std::cout << PX.col(NSTATE).transpose() << "\n";
        if (true) { // Normalize quaternions
          Eigen::Vector4d q = PX.block(6, NSTATE, 4, 1);
          PX.block(6, NSTATE, 4, 1) /= q.norm();
          q = xint.block(6, 0, 4, 1);
          xint.block(6, 0, 4, 1) /= q.norm();
        }
      }
      if (ii > 0) {
        // Perform kalman updates if needed:
        if ((int)ii % ngps == 0) {
          //std::cout << "GPS UPDATE!\n";
          KalmanGPSUpdate(xint, ngps * dt, &PX);
        }
        if ((int)ii % nimu == 0) {
          //std::cout << "IMU UPDATE!\n";
          KalmanIMUUpdate(xint, nimu * dt, &PX);
        }
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
                << qds_.lower_bound(t)->second.coeffs().transpose() << " "
                << PX.col(NSTATE).transpose() << "\n";
      //std::cout << t+dt << " " << PX.col(NSTATE).transpose() << "\n";
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
  void KalmanPredict(const Eigen::Array<double, NSTATE, NSTATE + 1> &PX,
                     Eigen::Array<double, NSTATE, NSTATE + 1> &PXdot,
                     double t) {
    MatrixNNd P = PX.block(0, 0, NSTATE, NSTATE);
    State x = PX.col(NSTATE);
    State xdot;
    Step(x, x, xdot, t);
    PXdot.col(NSTATE) = xdot;
    MatrixNNd A =
        math::jacobian<double, NSTATE, NSTATE>([this, x, t](State xjac) {
                                                 State xdot;
                                                 this->Step(xjac, x, xdot, t);
                                                 return xdot;
                                               },
                                               x);
    MatrixNNd Q = MatrixNNd::Zero();
    Q.diagonal() << 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    PXdot.block(0, 0, NSTATE, NSTATE) = A * P + P * A.transpose() + Q;
  }
  void KalmanGPSUpdate(const State &x, double dt,
                       Eigen::Array<double, NSTATE, NSTATE + 1> *PX) {
    if (PX == nullptr) {
      return;
    }
    // Actuall perform measurement
    Sensors::GPSMat zgps, zxhat;
    Sensors::GPSLin Cgps;
    Sensors::GPSCov Rgps;
    // Actually perform measurement.
    sensors_.GetGPS(x, dt, &zgps, nullptr, &Rgps);
    // Compute expected measurement
    sensors_.GetGPSNoiseless(PX->col(NSTATE), &zxhat, &Cgps);
    MatrixNNd Pmat = PX->block(0, 0, NSTATE, NSTATE);
    Eigen::Matrix<double, NSTATE, Sensors::NGPS> K =
        Pmat * Cgps.transpose() *
        (Cgps * Pmat * Cgps.transpose() + Rgps).inverse();
    PX->block(0, 0, NSTATE, NSTATE) = (Pmat.Identity() - K * Cgps) * Pmat;
    PX->col(NSTATE) = PX->col(NSTATE) + (K * (zgps - zxhat)).array();
  }
  void KalmanIMUUpdate(const State &x, double dt,
                       Eigen::Array<double, NSTATE, NSTATE + 1> *PX) {
    if (PX == nullptr) {
      return;
    }
    // Actuall perform measurement
    Sensors::IMUMat zimu, zxhat;
    Sensors::IMULin Cimu;
    Sensors::IMUCov Rimu;
    // Actually perform measurement.
    sensors_.GetIMU(x, dt, &zimu, nullptr, &Rimu);
    // Compute expected measurement
    sensors_.GetIMUNoiseless(PX->col(NSTATE), &zxhat, &Cimu);
    MatrixNNd Pmat = PX->block(0, 0, NSTATE, NSTATE);
    Eigen::Matrix<double, NSTATE, Sensors::NIMU> K =
        Pmat * Cimu.transpose() *
        (Cimu * Pmat * Cimu.transpose() + Rimu).inverse();
    PX->block(0, 0, NSTATE, NSTATE) = (Pmat.Identity() - K * Cimu) * Pmat;
    PX->col(NSTATE) = PX->col(NSTATE) + (K * (zimu - zxhat)).array();
  }
  void Step(const State &x, const State &xhat, State &xdot, double t) {
    Array3d pos_d, vel_d, accel_d, jerk_d;
    double psi_d, psidot_d;
    traj_.StateAtTime(t, &pos_d, &vel_d, &accel_d, &jerk_d, &psi_d, &psidot_d);

    Array3d moment_des;
    double thrust_des;
    Quaterniond qd_out;
    GoalForces(xhat, pos_d, vel_d, accel_d, jerk_d, psi_d, psidot_d, p_,
               &thrust_des, &moment_des, &qd_out);
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
