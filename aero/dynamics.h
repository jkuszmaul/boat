#pragma once

#include <Eigen/Core>

#include "params.h"
#include <random>

namespace aero {

constexpr int NSTATE = 13;
typedef Eigen::Array<double, NSTATE, 1> State;
typedef Eigen::Array<double, NSTATE, NSTATE> ArrayNNd;
typedef Eigen::Matrix<double, NSTATE, NSTATE> MatrixNNd;
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

class Sensors {
 public:
  static constexpr int NGPS = 6, NIMU = 10;
  // GPSMat will contain (pnorth, peast, pdown, vnorth, veast, vdown)
  typedef Eigen::Matrix<double, NGPS, 1> GPSMat;
  typedef Eigen::Matrix<double, NGPS, NGPS> GPSCov;
  typedef Eigen::Matrix<double, NGPS, NSTATE> GPSLin;
  // IMUMat will contain (accelx, accely, accelz, gyrox, gyroy, gyroz, magx, magy, magz, pressurez)
  // Accel in m/s^2
  // Gyro in rad/s
  // Magnetometer returns direction of true north relative to us
  // Pressure returns pressure difference between current position and origin
  typedef Eigen::Matrix<double, NIMU, 1> IMUMat;
  typedef Eigen::Matrix<double, NIMU, NIMU> IMUCov;
  typedef Eigen::Matrix<double, NIMU, NSTATE> IMULin;
  Sensors();
  void GetGPSNoiseless(State x, GPSMat *gps, GPSLin *Cgps);
  void GetIMUNoiseless(State x, IMUMat *imu, IMULin *Cimu);
  void GetGPS(State x, double dt, GPSMat *gps, GPSLin *Cgps, GPSCov *Rgps);
  void GetIMU(State x, double dt, IMUMat *imu, IMULin *Cimu, IMUCov *Rimu);

 private:
  IMUMat hIMU(State x);

  template <int N>
  Eigen::Matrix<double, N, 1>  GaussianVec() {
    Eigen::Matrix<double, N, 1> x;
    for (int ii = 0; ii < N; ++ii) {
      x(ii, 0) = normal_(mersenne_);
    }
    return x;
  }
  // Model parameters:
  Eigen::Vector3d kgps{1. / 1100., 1. / 1100., 1. / 1100.},          // s^{-1}
      stdgpspos{0.21, 0.21, 0.4},                                    // m
      stdgpsvel{0.01, 0.01, 0.05},                                   // m
      stdmag{0.3 * M_PI / 180., .3 * M_PI / 180., .3 * M_PI / 180.}, // rad
      stdgyro{0.13 * M_PI / 180, 0.13 * M_PI / 180,
              0.13 * M_PI / 180},                        // rad / s
      stdaccel{0.0025 * 9.8, 0.0025 * 9.8, 0.0025 * 9.8} // m / s^2
  ;
  double rho_ = 1.225; // kg / m^3
  double stdpres{0.01}; // kPa
  // Last noise from GPS
  GPSMat gps_noise_;
  // Last noise from IMU
  IMUMat imu_noise_;

  // For making random noise:
  std::random_device rd_;
  std::mt19937 mersenne_{rd_()};
  std::normal_distribution<> normal_;
};  // class Sensors

}  // namespace aero
