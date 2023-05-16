// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LOGGER_H
#define UNO_LOGGER_H

#include <iostream>

#define RED    "\x1B[31m"
// #define GREEN   "\x1B[32m"
#define YELLOW  "\x1B[33m"
// #define BLUE    "\x1B[34m"
// #define MAGENTA "\x1B[35m"
// #define CYAN    "\x1B[36m"
// #define WHITE   "\x1B[37m"
#define RESET "\x1B[0m"

enum Level {
    ERROR = 0, WARNING, INFO, DEBUG, DEBUG2
};

class Logger {
public:
    static Level level;
    static void set_logger(const std::string& logger_level);
};

template <typename T>
const Level& operator<<(const Level& level, T& element) {
    if (level <= Logger::level) {
        std::cout << element;
    }
    return level;
}

template <typename T>
const Level& operator<<(const Level& level, const T& element) {
    if (level <= Logger::level) {
        std::cout << element;
    }
    return level;
}

inline void Logger::set_logger(const std::string& logger_level) {
   if (logger_level == "ERROR") {
      Logger::level = ERROR;
   }
   else if (logger_level == "WARNING") {
      Logger::level = WARNING;
   }
   else if (logger_level == "INFO") {
      Logger::level = INFO;
   }
   else if (logger_level == "DEBUG") {
      Logger::level = DEBUG;
   }
   else if (logger_level == "DEBUG2") {
      Logger::level = DEBUG2;
   }
   else {
      throw std::out_of_range("The logger level " + logger_level + " was not found");
   }
}

#endif // UNO_LOGGER_H