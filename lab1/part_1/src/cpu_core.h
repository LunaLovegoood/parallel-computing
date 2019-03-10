#ifndef CPU_CORE_H_
#define CPU_CORE_H_

#include "process.h"

#include <mutex>

class cpu_core {
public:
  bool execute(const process_t &process) const;

private:
  mutable bool is_busy_{ false };
  mutable std::mutex mutex_;
};

#endif
