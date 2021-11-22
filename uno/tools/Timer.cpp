#include "Timer.hpp"
#include <chrono>
#include <ctime>

void Timer::start() {
   this->start_time = std::clock();
}

void Timer::stop() {
   this->end_time = std::clock();
}

double Timer::get_duration() const {
   return static_cast<double>(this->end_time - this->start_time) / static_cast<double>(CLOCKS_PER_SEC);
}

char* Timer::get_current_date() {
   const auto current_time = std::chrono::system_clock::now();
   const auto current_time_2 = std::chrono::system_clock::to_time_t(current_time);
   return std::ctime(&current_time_2);
}