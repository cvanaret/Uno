// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STRINGS_H
#define UNO_STRINGS_H

#include <string_view>

namespace uno {
   // case-insensitive comparison of ASCII strings
   bool case_insensitive_equals(std::string_view a, std::string_view b);
} // namespace

#endif // UNO_STRINGS_H