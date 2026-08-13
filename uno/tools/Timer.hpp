// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TIMER_H
#define UNO_TIMER_H

#include <chrono>
#include <string>

namespace uno {
   class Timer {
   public:
      Timer();
      [[nodiscard]] double get_duration() const;
      [[nodiscard]] static std::string get_current_date();

   private:
      std::chrono::time_point<std::chrono::steady_clock> start;
   };
} // namespace

#endif //UNO_TIMER_H