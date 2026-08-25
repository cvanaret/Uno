// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Timer.hpp"
#include <chrono>
#include <iomanip>
#include <sstream>

namespace uno {
   void Timer::start() {
      this->start_point = std::chrono::steady_clock::now();
   }

   void Timer::stop() {
      const auto now = std::chrono::steady_clock::now();
      this->duration += std::chrono::duration<double>(now - this->start_point).count();
   }

   double Timer::get_elapsed_time() const {
      const auto now = std::chrono::steady_clock::now();
      return std::chrono::duration<double>(now - this->start_point).count();
   }

   double Timer::get_duration() const {
      return this->duration;
   }

   std::string Timer::get_current_date() {
      const auto now = std::chrono::system_clock::now();
      const std::time_t time = std::chrono::system_clock::to_time_t(now);

      std::tm calendar_time{};
#if defined(_WIN32)
      localtime_s(&calendar_time, &time);       // note: reversed argument order vs POSIX
#else
      localtime_r(&time, &calendar_time);
#endif
      std::ostringstream stream;
      stream << std::put_time(&calendar_time, "%a %b %d %H:%M:%S %Y");
      return stream.str();
   }
} // namespace