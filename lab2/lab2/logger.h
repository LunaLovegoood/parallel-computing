#ifndef LOGGER_H_
#define LOGGER_H_

#include <mutex>
#include <string>
#include <fstream>

class Logger {
 public:
  static Logger& getInstance();

  void log(const std::string &msg);

  Logger(const Logger &logger) = delete;
  Logger(Logger &&logger) = delete;

  Logger& operator=(const Logger& logger) = delete;
  Logger& operator=(Logger&& logger) = delete;

 private:
  Logger();
  std::mutex mutex_;

  std::ofstream out_;
};

#endif  // LOGGER_H_
