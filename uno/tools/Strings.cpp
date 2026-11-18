// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cctype>
#include "Strings.hpp"

namespace uno {
   bool case_insensitive_equals(std::string_view a, std::string_view b) {
      return a.size() == b.size() &&
         std::equal(a.begin(), a.end(), b.begin(), [](unsigned char x, unsigned char y) {
            return std::toupper(x) == std::toupper(y);
         });
   }

   std::string abbreviate(std::string_view string) {
      constexpr size_t max_printed_length = 64;
      if (string.size() <= max_printed_length) {
         return std::string(string);
      }
      // do not cut a UTF-8 character: skip back over its continuation bytes (10xxxxxx)
      size_t length = max_printed_length;
      while (0 < length && (static_cast<unsigned char>(string[length]) & 0xC0) == 0x80) {
         --length;
      }
      return std::string(string.substr(0, length)) + "... (" + std::to_string(string.size()) + " bytes)";
   }
} // namespace
