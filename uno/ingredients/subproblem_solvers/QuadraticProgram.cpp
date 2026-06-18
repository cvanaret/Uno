// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QuadraticProgram.hpp"

namespace uno {
   QuadraticProgram::QuadraticProgram(size_t number_variables, size_t number_constraints):
      number_variables(number_variables), number_constraints(number_constraints) {
   }

   // out-of-line anchor for the vtable
   QuadraticProgram::~QuadraticProgram() = default;
} // namespace
