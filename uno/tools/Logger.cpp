// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Logger.hpp"

namespace uno {
   void Logger::set_logger(const std::string& logger_level) {
      if (logger_level == "SILENT") {
         Logger::level = SILENT;
      }
      else if (logger_level == "DISCRETE") {
         Logger::level = DISCRETE;
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
      else if (logger_level == "DEBUG3") {
         Logger::level = DEBUG3;
      }
      else {
         throw std::out_of_range("The logger level " + logger_level + " was not found");
      }
   }
} // namespace