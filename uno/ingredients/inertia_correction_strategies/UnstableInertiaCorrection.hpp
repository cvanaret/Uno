// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNSTABLEINERTIACORRECTION_H
#define UNO_UNSTABLEINERTIACORRECTION_H

#include <exception>

namespace uno {
   struct UnstableInertiaCorrection: public std::exception {

      [[nodiscard]] const char* what() const noexcept override {
         return "The inertia correction got unstable (delta_w > threshold)";
      }
   };
} // namespace

#endif // UNO_UNSTABLEINERTIACORRECTION_H