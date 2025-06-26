// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LITERALS_H
#define UNO_LITERALS_H

#include <array>

namespace uno {
   // make_literal_array creates an array of literals (const char*)
   // the return type is correct, even if there is no argument supplied
   template <typename... Types>
   constexpr std::array<const char*, sizeof...(Types)> make_literal_array(Types&&... types) {
      return {types...};
   }
} // namespace

#endif // UNO_LITERALS_H