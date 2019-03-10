#include "cpu_core.h"

#include <thread>

bool cpu_core::execute(const process_t &process) const {
  auto execution_function = [this](const process_t &process) {
    process();
    std::lock_guard<std::mutex> lock(mutex_);
    is_busy_ = false;
  };

  std::lock_guard<std::mutex> lock(mutex_);
  if (!is_busy_) {
    is_busy_ = true;
    std::thread execution_thread(execution_function, process);
    execution_thread.detach();
    return true;
  }

  return false;
}
