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
} // namespace