#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <chrono>
#include <string>

std::mutex first_m;
std::mutex second_m;

std::condition_variable first_var;
std::condition_variable second_var;

bool first_notified = false;
bool second_notified = false;

void print_three(std::string first, std::string second, std::string third);

int main() {
  auto first = std::thread(
    [&](std::string first, std::string second, std::string third) {
      std::unique_lock<std::mutex> lock(first_m);
      print_three(first, second, third);
      first_notified = true;
      first_var.notify_one();
    },
    "", "Khomiak", ""
  );

  auto second = std::thread(
    [&](std::string first, std::string second, std::string third) {
      std::unique_lock<std::mutex> first_lock(first_m);
      while (!first_notified) {
        first_var.wait(first_lock);
      }

      std::unique_lock<std::mutex> second_lock(second_m);
      print_three(first, second, third);
      second_notified = true;
      second_var.notify_one();
    },
    "", "Yurii", ""
  );

  auto third = std::thread(
    [&](std::string first, std::string second, std::string third) {
      std::unique_lock<std::mutex> lock(second_m);
      while (!second_notified) {
        second_var.wait(lock);
      }
      
      print_three(first, second, third);
    },
    "", "Vasyliovych", ""
  );

  first.join();
  second.join();
  third.join();

  return 0;
}

void print_three(std::string first, std::string second, std::string third) {
  std::cout << first << ' ' << second << ' ' << third << '\n';
}
