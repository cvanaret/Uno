// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Parameterization.hpp"

namespace uno {
   void Parameterization::set(std::string_view key, double value) {
      this->parameters[key] = value;
   }

   double Parameterization::get(std::string_view key) const {
      return this->parameters.at(key);
   }
} // namespace