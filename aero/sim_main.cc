#include "sim.h"
#include "gflags.h"
#include "glog/logging.h"

int main(int argc, char *argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  aero::Simulate sim;
  sim.Run(7.25);
  //sim.Run(0.2);
}
