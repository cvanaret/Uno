// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERSPARSEREPRESENTATION_H
#define UNO_LINEARSOLVERSPARSEREPRESENTATION_H

#include "SolverWorkspace.hpp"
#include "../interfaces/C/uno_int.h"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class LinearSolverSparseRepresentation: public SolverWorkspace {
   public:
      LinearSolverSparseRepresentation() = default;
      ~LinearSolverSparseRepresentation() override = default;

      // TODO: these shouldn't be here!
      std::vector<uno_int> jacobian_row_indices{};
      std::vector<uno_int> jacobian_column_indices{};

      size_t dimension{0};
      size_t number_nonzeros{0};
      Vector<double> matrix_values{};
      Vector<double> rhs{};
      Vector<double> solution{};

      virtual void initialize_hessian(const Subproblem& subproblem) = 0;
      virtual void initialize_augmented_system(const Subproblem& subproblem) = 0;
   };
} // namespace

#endif // UNO_LINEARSOLVERSPARSEREPRESENTATION_H