// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Timer.hpp"
#include <chrono>
#include <ctime>

namespace uno {
   Timer::Timer(): start_time(std::clock()) {
   }

   double Timer::get_duration() const {
      return static_cast<double>(std::clock() - this->start_time) / static_cast<double>(CLOCKS_PER_SEC);
   }

   char* Timer::get_current_date() {
      const auto current_time = std::chrono::system_clock::now();
      const auto formatted_current_time = std::chrono::system_clock::to_time_t(current_time);
      return std::ctime(&formatted_current_time);
   }
} // namespace
