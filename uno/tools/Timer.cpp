// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

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
   const auto formatted_current_time = std::chrono::system_clock::to_time_t(current_time);
   return std::ctime(&formatted_current_time);
}