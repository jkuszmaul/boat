#include "dynamics.h"

#include "math/jacobian.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

DEFINE_double(kp, 4.0, "kP");
DEFINE_double(kd, 2.0, "kD");
DEFINE_double(kpm, 10.0, "kPM");
DEFINE_double(kdm, 0.1, "kDM");

namespace aero {

using Eigen::Quaterniond;

Sensors::Sensors() {}

void Sensors::GetGPSNoiseless(State x, GPSMat *gps, GPSLin *Cgps) {
  if (gps != nullptr) {
    *gps = x.block(0, 0, 6, 1);
  }
  if (Cgps != nullptr) {
    *Cgps = GPSLin::Identity();
  }
}

void Sensors::GetIMUNoiseless(State x, IMUMat *imu, IMULin *Cimu) {
  if (imu != nullptr) {
    *imu = hIMU(x);
  }
  if (Cimu != nullptr) {
    *Cimu = math::jacobian<double, NSTATE, NIMU>(
        std::bind(&Sensors::hIMU, this, std::placeholders::_1), x);
  }
}

void Sensors::GetGPS(State x, double dt, Sensors::GPSMat *gps,
                     Sensors::GPSLin *Cgps, Sensors::GPSCov *Rgps) {
  GetGPSNoiseless(x, gps, Cgps);

  // Now, calculate noise:
  if (Rgps != nullptr) {
    // Compute covariance matrix as being the sum of semi-biased and white noise
    // characteristics:
    *Rgps = Rgps->Zero();
    // Values pulled from book. Not sure how accurate velocity is.
    Rgps->diagonal() << 6.6, 6.6, 9.2, 0.05, 0.05, 0.1;
    *Rgps = Rgps->cwiseAbs2(); // Square for covariance
  }


  if (gps != nullptr) {
    Eigen::Vector3d randpos = GaussianVec<3>(), randvel = GaussianVec<3>();
    // Compute actual noise, from book:
    gps_noise_.block(0, 0, 3, 1) =
        Eigen::exp(-kgps.array() * dt).matrix().cwiseProduct(
            gps_noise_.block(0, 0, 3, 1)) +
        stdgpspos.cwiseProduct(randpos);
    gps_noise_.block(3, 0, 3, 1) = stdgpsvel.cwiseProduct(randvel);
    *gps += gps_noise_;
  }
}

void Sensors::GetIMU(State x, double dt, Sensors::IMUMat *imu,
                     Sensors::IMULin *Cimu, Sensors::IMUCov *Rimu) {
  GetIMUNoiseless(x, imu, Cimu);
  Eigen::Matrix<double, 10, 1> stds;
  stds << stdaccel, stdgyro, stdmag, stdpres;
  if (Rimu != nullptr) {
    *Rimu = Rimu->Zero();
    Rimu->diagonal() = stds.cwiseAbs2();
  }

  if (imu != nullptr) {
    Eigen::Matrix<double, 10, 1> randn = GaussianVec<10>();
    *imu += randn.cwiseProduct(stds);
  }
}

Sensors::IMUMat Sensors::hIMU(State x) {
  IMUMat imu = IMUMat::Zero();
  Quaterniond q;
  q.coeffs() = x.block(6, 0, 4, 1);
  Eigen::Matrix3d R = q.normalized().toRotationMatrix(); // Body->Inertial
  Eigen::Matrix3d Rinv = R.transpose(); // Inertial->Body
  Eigen::Vector3d Omega = x.block(10, 0, 3, 1);
  Eigen::Vector3d vel = x.block(3, 0, 3, 1);
  // Accelerometer: Centripetal accel + gravity:
  imu.block(0, 0, 3, 1) =
      Omega.cross(Rinv * vel) + Rinv * Eigen::Vector3d(0, 0, 1);
  // Gyro:
  imu.block(3, 0, 3, 1) = Omega;
  // Magnetometer, model as returning the current apparent direction of
  // true north
  imu.block(6, 0, 3, 1) = Rinv * Eigen::Vector3d(1, 0, 0);
  // Barometer, presssure relative to surface
  // Given relatively small altitude differences, use constant air density assumption
  imu(10, 0) = -rho_ * 9.8 * x(2, 0);
  return imu;
}

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
