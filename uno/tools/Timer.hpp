// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TIMER_H
#define UNO_TIMER_H

#include <ctime>

namespace uno {
   // timer starts upon creation
   class Timer {
   public:
      Timer();
      [[nodiscard]] double get_duration() const;
      [[nodiscard]] static char* get_current_date();

   private:
      std::clock_t start_time;
   };
} // namespace

#endif //UNO_TIMER_H
