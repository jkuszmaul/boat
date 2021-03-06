#pragma once
#include <map>
#include <condition_variable>
#include <shared_mutex>
#include <mutex>
#include <atomic>
#include <string>

namespace sailbot {
namespace testing {

class TestQueue {
 public:
  TestQueue(const std::string &name, size_t max_size)
      : name_(name), max_size_(max_size) {
    std::unique_lock<std::mutex> lck_(map_lock_);
    conditions_[name_].cnt++;
  }
  ~TestQueue() {
    std::unique_lock<std::mutex> lck_(map_lock_);
    conditions_[name_].cnt--;
    if (conditions_[name_].cnt.load() == 0) {
      conditions_.erase(name_);
    }
  }
  void send(const void *msg, size_t size);
  bool receive(void *msg, size_t size, size_t &rcvd);
 private:
  struct QueueData {
    std::condition_variable cond;
    std::unique_ptr<uint8_t[]> data;
    std::atomic<int> cnt{0};
    size_t data_len = 0;
    int inc = 0;
  };
  static std::mutex map_lock_;
  static std::map<std::string, QueueData> conditions_;

  std::string name_;
  const size_t max_size_;
  // Assumes that last_inc_ is called in only one thread.
  int last_inc_ = 0;
};

}  // testing
}  // sailbot
