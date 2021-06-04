#include "Timer.hpp"

void Timer::start() {
   this->start_time = std::clock();
}

void Timer::stop() {
   this->end_time = std::clock();
}

double Timer::get_time() const {
   return (double) (this->end_time - this->start_time) / (double) CLOCKS_PER_SEC;
}