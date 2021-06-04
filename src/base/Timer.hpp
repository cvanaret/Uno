#ifndef TIMER_H
#define TIMER_H

#include <ctime>

class Timer {
public:
   Timer() = default;
   void start();
   void stop();
   [[nodiscard]] double get_time() const;

private:
   std::clock_t start_time, end_time;
};

#endif //TIMER_H
