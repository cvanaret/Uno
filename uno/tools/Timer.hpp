// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TIMER_H
#define UNO_TIMER_H

#include <chrono>
#include <string>

namespace uno {
   // timer that accumulates a total duration (gaps between start() and stop())
   class Timer {
   public:
      Timer() = default;

      void start();
      void stop();
      [[nodiscard]] double get_elapsed_time() const;
      [[nodiscard]] double get_duration() const;

      [[nodiscard]] static std::string get_current_date();

   private:
      std::chrono::steady_clock::time_point start_point;
      double duration{0.};
   };
} // namespace

#endif //UNO_TIMER_H