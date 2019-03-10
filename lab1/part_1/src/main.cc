#include "cpu_core.h"
#include "process.h"

#include <iostream>
#include <array>
#include <chrono>
#include <thread>
#include <cstdint>

constexpr uint64_t NUMBER_OF_PROCESSES = 1000;
constexpr uint64_t NUMBER_OF_CPUS = 2;

struct run_stats {
  uint64_t successful = 0;
  uint64_t destroyed = 0;
};

void simulate_process_flow();
bool execute_process(
  const std::array<cpu_core, NUMBER_OF_CPUS> &cpu_cores,
  const process_t &process
);

void update_counters(run_stats *stats, bool isDestroyed);
void delay(double delay_in_ms);

std::ostream& operator<<(std::ostream& stream, const run_stats &stats);

int main() {
  std::string(foo);
  simulate_process_flow();
  return 0;
}

void simulate_process_flow() {
  std::array<cpu_core, NUMBER_OF_CPUS> cpu_cores{};
  run_stats stats{};
  
  for (int i = 0; i < NUMBER_OF_PROCESSES; i++) {
    process_t process;

    bool is_executed = execute_process(cpu_cores, process);
    update_counters(&stats, is_executed);

    if (stats.destroyed % 50) {
      delay(0.9);
    }
  }

  std::cout << stats;
}

bool execute_process(
  const std::array<cpu_core, NUMBER_OF_CPUS> &cpu_cores,
  const process_t &process
) {
  for (int i = 0; i < cpu_cores.size(); i++) {
    if (cpu_cores[i].execute(process)) {
      return true;
    }
  }
  return false;
}

void update_counters(run_stats *stats, bool is_executed) {
  if (is_executed) {
    stats->successful++;
  } else {
    stats->destroyed++;
  }
}

void delay(double delay_in_ms) {
  std::this_thread::sleep_for(
    std::chrono::duration<double, std::milli>(delay_in_ms)
  );
}

std::ostream& operator<<(std::ostream& stream, const run_stats &stats) {
  std::cout << "Successful: " << stats.successful << '\n';
  std::cout << "Destroyed: " << stats.destroyed << '\n';
  return stream;
}
