#pragma once
#include "util/node.h"
#include "control/actuator_cmd.pb.h"
#include <mutex>

namespace sailbot {
namespace control {

class SimpleControl : public Node {
 public:
  SimpleControl();

  void Iterate() override;
 private:
  constexpr static float dt = 0.01;
  msg::SailCmd *sail_msg_;
  msg::RudderCmd *rudder_msg_;
  msg::BoatState *boat_state_;
  std::mutex boat_state_mutex_;
  std::atomic<float> wind_x_{0}, wind_y_{0};
  ProtoQueue<msg::SailCmd> sail_cmd_;
  ProtoQueue<msg::RudderCmd> rudder_cmd_;

  // Various useful bits
  float last_goal_ = 0;
  float goal_cost_ = 0;
};

}  // control
}  // sailbot