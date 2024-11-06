// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LOGGER_H
#define UNO_LOGGER_H

#include <string>
#include <iostream>

namespace uno {
   enum Level {
       SILENT = 0, DISCRETE, WARNING, INFO, DEBUG, DEBUG2, DEBUG3
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
} // namespace

#endif // UNO_LOGGER_H
