#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "dynamics.h"
#include <vector>

namespace aero {
void GoalForces(const State &x, const Array3d &r_d, const Array3d &rdot_d,
                const Array3d &rddot_d, const Array3d &rdddot_d,
                const double psi_d, const double psidot_d, const Param &p,
                double *thrust_des, Array3d *moment_des,
                Eigen::Quaterniond *qd_out = nullptr);

class Trajectory {
 public:
  Trajectory(std::vector<double> a, std::vector<double> b,
             std::vector<double> c)
       : a_(a), b_(b), c_(c) {}

  void StateAtTime(double t, Array3d *pos, Array3d *vel, Array3d *accel,
                   Array3d *jerk, double *psi, double *psidot) {
    if (pos != nullptr) {
      PolyPos(t, 0, pos);
    }
    if (vel != nullptr) {
      PolyPos(t, 1, vel);
    }
    if (accel != nullptr) {
      PolyPos(t, 2, accel);
    }
    if (jerk != nullptr) {
      PolyPos(t, 3, jerk);
    }
    if (psi != nullptr) {
      double pydot = Poly(t, b_, 1);
      double pxdot = Poly(t, a_, 1);
      *psi = -std::atan2(pydot, pxdot);
    }
    if (psidot != nullptr) {
      double pydot = Poly(t, b_, 1);
      double pxdot = Poly(t, a_, 1);
      double pyddot = Poly(t, b_, 2);
      double pxddot = Poly(t, a_, 2);
      double vel2 = pxdot * pxdot + pydot * pydot;
      if (vel2 > 1e-4) {
        *psidot =
            -(pxdot * pyddot - pxddot * pydot) / vel2;
      } else {
        *psidot = 0.0;
      }
    }
  }
  void PolyPos(double t, int deriv, Array3d *p) {
    *p << Poly(t, a_, deriv), Poly(t, b_, deriv), Poly(t, c_, deriv);
  }
 private:
  // Compute the deriv'th derivative of a polynomial with coefficients a at x
  static double Poly(double x, std::vector<double> a, int deriv) {
    double retval = 0;
    double accum = 1;
    for (int ii = deriv; ii < a.size(); ++ii) {
      double mult = 1;
      for (int jj = 0; jj < deriv; ++jj) {
        mult *= (ii - jj);
      }
      retval += mult * accum * a[ii];
      accum *= x;
    }
    return retval;
  }
  std::vector<double> a_, b_, c_;
};  // class Trajectory
}  // namespace aero
