#include "Timer.hpp"

void Timer::start() {
   this->start_time = std::clock();
}

void Timer::stop() {
   this->end_time = std::clock();
}

double Timer::get_time() const {
   return static_cast<double>(this->end_time - this->start_time) / static_cast<double>(CLOCKS_PER_SEC);
}