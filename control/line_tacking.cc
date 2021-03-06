#include "line_tacking.h"
#include "control/actuator_cmd.pb.h"
#include "util.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace sailbot {
namespace control {

constexpr int LineTacker::N_WAYPOINTS;

LineTacker::LineTacker()
    : Node(0.01), kCloseHaul(M_PI / 4), kWindTol(0.1),
      cur_pos_({0, 0}), wind_dir_(0), cur_theta_(0),
      heading_msg_(AllocateMessage<msg::HeadingCmd>()),
      heading_cmd_("heading_cmd", true) {
  RegisterHandler<msg::Vector3f>("wind", [this](const msg::Vector3f &msg) {
    wind_dir_ = std::atan2(msg.y(), msg.x());
  });
  // TODO(james): Thread-safeness.
  RegisterHandler<msg::BoatState>("boat_state",
                                  [this](const msg::BoatState &msg) {
    Eigen::Quaternionf orientation(msg.orientation().w(), msg.orientation().x(),
                                   msg.orientation().y(),
                                   msg.orientation().z());
    Eigen::Matrix3f rot = orientation.toRotationMatrix();
    cur_theta_ = std::atan2(rot(1, 0), rot(0, 0));
    cur_pos_.x = msg.pos().x();
    cur_pos_.y = msg.pos().y();
  });

  RegisterHandler<msg::WaypointList>("waypoints",
                                     [this](const msg::WaypointList &msg) {
    // TODO(james): Thread-safety
    int starti = msg.restart() ? 0 : msg_i_offset_ + i_;
    int cnt = std::min(msg.points_size()-starti, N_WAYPOINTS);
    for (int i = 0; i < cnt; ++i) {
      waypoints_[i] = {msg.points(i + starti).x(), msg.points(i + starti).y()};
    }
    i_ = 0;
    msg_i_offset_ = starti;
    way_len_ = cnt;
  });
}

void LineTacker::Iterate() {
  if (way_len_ == 0) {
    LOG(INFO) << "No Waypoints--not doing anything";
    return;
  }
  double gh = GoalHeading();
  heading_msg_->set_heading(gh);
  heading_cmd_.send(heading_msg_);
}

/**
 * Basic algorithm:
 * If we have a straight-line path from our current position to the goal
 * position and the straight line isn't too close to being close-hauled, then
 * aim to go straight there.
 * Otherwise, depends on the state:
 *   -If not over/at the line boundary, go to whichever closehaul is closest.
 *   -If over/at the boundary, to closehaul that brings us into the line.
 */
float LineTacker::GoalHeading() {
  Point start = waypoints_[i_];
  Point end = waypoints_[i_+1];
  float upwind = util::norm_angle(wind_dir_ - M_PI);
  float min_closehaul = util::norm_angle(upwind - kCloseHaul);
  float max_closehaul = util::norm_angle(upwind + kCloseHaul);
  float dy = end.y - cur_pos_.y;
  float dx = end.x - cur_pos_.x;
  float nominal_heading = std::atan2(dy, dx);
  float dist = std::sqrt(dy * dy + dx * dx);
  float goal_wind_diff = util::norm_angle(nominal_heading - upwind);
  if (dist < 2 && i_ < way_len_ - 2) {
    ++i_;
  }

  if (std::abs(goal_wind_diff) > kCloseHaul + kWindTol) {
    // It is worth it to go ahead and go straight there.
    return nominal_heading;
  }

  float dist_to_line = DistanceFromLine(start, end, cur_pos_);

  if (std::abs(dist_to_line) >= bounds_[i_]) {
    // Too far from the path, so go back.
    if (util::norm_angle(wind_dir_ - nominal_heading) > 0) {
      return max_closehaul;
    } else {
      return min_closehaul;
    }
  } else {
    // Close enough to the line, so go wherever is closest.
    if (ApparentWind() > 0) {
      return max_closehaul;
    } else {
      return min_closehaul;
    }
  }

  return nominal_heading;
}

float LineTacker::DistanceFromLine(Point start, Point end, Point loc) {
  // Account for instances when we are near singularities:
  if (std::abs(end.y - start.y) < 1e-3) {
    return (start.y + end.y) / 2 - loc.y;
  } else if (std::abs(end.x - start.x) < 1e-3) {
    return (start.x + end.x) / 2 - loc.x;
  }
  // Draw a line that goes through loc and is perpendicular to the original.
  // Express it as y = perp_slope * x + perp_intercept;
  float perp_slope = -(end.x - start.x) / (end.y - start.y);
  float perp_intercept = loc.y - perp_slope * loc.x;
  float orig_slope = -1 / perp_slope;
  float orig_intercept = end.y - orig_slope * end.x;
  // Now, find intersection between perp_line and original line.
  float x = (perp_intercept - orig_intercept) / (orig_slope - perp_slope);
  float y = perp_slope * x + perp_intercept;
  float dx = x - loc.x;
  float dy = y - loc.y;
  float dist = std::sqrt(dx * dx + dy * dy);
  dist *= dx > 0 ? 1 : -1;
  return dist;
}

float LineTacker::ApparentWind() {
  return util::norm_angle(wind_dir_ - cur_theta_);
}

}  // namespace control
}  // namespace sailbot
