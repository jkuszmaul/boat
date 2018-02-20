#include "eigen_int.hpp"
#include <iostream>

using namespace boost::numeric::odeint;
using namespace Eigen;

typedef Array<double, 3, 1> Arrayd;
void f(const Arrayd &x, Arrayd &xdot, const double t) {
  xdot = -Eigen::log(Eigen::abs(x));
}

int main(int argc, char *argv[]) {
  dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<Arrayd>>>
      integrator;
  const Arrayd x0(0.01, 0.001, 0.01);
  integrator.initialize(x0, 0, 10000);
  std::pair<double, double> tstep = integrator.do_step(&f);
  std::cout << "t0: " << tstep.first << " tf: " << tstep.second << std::endl;
  Arrayd x;
  for (double ii = 0; ii < 10; ii += 0.01) {
    integrator.calc_state(ii, x);
    std::cout << ii << ": " << x.transpose() << std::endl;
  }
}
