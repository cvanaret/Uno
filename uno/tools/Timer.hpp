// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_TIMER_H
#define UNO_TIMER_H

#include <ctime>

class Timer {
public:
   Timer() = default;
   void start();
   void stop();
   [[nodiscard]] double get_duration() const;
   [[nodiscard]] static char* get_current_date();

private:
   std::clock_t start_time, end_time;
};

#endif //UNO_TIMER_H
