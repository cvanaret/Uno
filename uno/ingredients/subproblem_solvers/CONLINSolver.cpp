// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "CONLINSolver.hpp"
#include "optimization/EvaluationCache.hpp"
#include <cmath>

namespace uno {
   CONLINSolver::CONLINSolver(size_t number_variables, size_t number_constraints, const Options& options) :
         SCPSolver(number_variables, number_constraints, options) {
   }

   void CONLINSolver::compute_diagonal_hessian(const Vector<double>& initial_point, Evaluations& current_evaluations, std::vector<double>& Dx) {
      // CONLIN Approximation
      // If df/dx > 0 -> Linear (D_x = 0)
      // If df/dx < 0 -> Reciprocal (D_x = -2 * df/dx / x_0)
      
      for (size_t j = 0; j < this->number_variables; ++j) {
         double df0dx = current_evaluations.objective_gradient[j];
         if (df0dx < 0.0) {
             Dx[j] -= 2.0 * df0dx / std::max(1e-8, initial_point[j]);
         }
      }

      size_t num_jacobian_nonzeros = current_evaluations.jacobian_values.size();
      for (size_t k = 0; k < num_jacobian_nonzeros; ++k) {
         size_t i = current_evaluations.jacobian_sparsity->row_indices[k];
         size_t j = current_evaluations.jacobian_sparsity->column_indices[k];
         if (j < this->number_variables) {
            double dfdx = current_evaluations.jacobian_values[k];
            if (dfdx < 0.0) {
               // Approximate multiplier weighting using constraints
               double d2f = -2.0 * dfdx / std::max(1e-8, initial_point[j]);
               Dx[j] += std::max(0.0, current_evaluations.constraints[i]) * d2f;
            }
         }
      }
   }
} // namespace
