// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STRINGS_H
#define UNO_STRINGS_H

#include <string>
#include <string_view>

namespace uno {
   // case-insensitive comparison of ASCII strings
   bool case_insensitive_equals(std::string_view a, std::string_view b);

   // a string that comes from the user (e.g. the name of an unknown option) in an error message: it may be arbitrarily
   // long, so only its first 64 bytes are shown (without cutting a UTF-8 character in half), followed by its length in bytes
   std::string abbreviate(std::string_view string);
} // namespace

#endif // UNO_STRINGS_H