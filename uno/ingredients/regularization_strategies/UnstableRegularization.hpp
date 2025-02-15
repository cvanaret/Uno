// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNSTABLEREGULARIZATION_H
#define UNO_UNSTABLEREGULARIZATION_H

#include <exception>

namespace uno {
   struct UnstableRegularization : public std::exception {

      [[nodiscard]] const char* what() const noexcept override {
         return "The inertia correction got unstable (delta_w > threshold)";
      }
   };
} // namespace

#endif // UNO_UNSTABLEREGULARIZATION_H