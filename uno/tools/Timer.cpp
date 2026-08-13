// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Timer.hpp"
#include <chrono>
#include <ctime>

namespace uno {
   // timer starts upon creation
   Timer::Timer(): start(std::chrono::steady_clock::now()) {
   }

   double Timer::get_duration() const {
      const auto now = std::chrono::steady_clock::now();
      return std::chrono::duration<double>(now - this->start).count();
   }

   char* Timer::get_current_date() {
      const auto current_time = std::chrono::system_clock::now();
      const auto formatted_current_time = std::chrono::system_clock::to_time_t(current_time);
      return std::ctime(&formatted_current_time);
   }
} // namespace