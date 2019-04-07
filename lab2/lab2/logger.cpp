#include "logger.h"

Logger::Logger() {
  out_.open("..\\out\\out.txt", std::ios::out);
}

Logger& Logger::getInstance() {
  static Logger instance;
  return instance;
}

void Logger::log(const std::string &msg) {
  std::lock_guard<std::mutex> lock(mutex_);
  out_ << msg << std::endl;
}
