// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Preprocessing.hpp"
#include "tools/Range.hpp"

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
void Preprocessing::compute_least_square_multipliers(const Model& model, SymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver& solver,
      Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_norm) {
   const size_t number_variables = current_iterate.primals.size();
   current_iterate.evaluate_objective_gradient(model);
   current_iterate.evaluate_constraint_jacobian(model);

   /******************************/
   /* build the symmetric matrix */
   /******************************/
   matrix.reset();
   // identity block
   for (size_t i = 0; i < number_variables; i++) {
      matrix.insert(1., i, i);
      matrix.finalize(i);
   }
   // Jacobian of general constraints
   for (size_t j = 0; j < model.number_constraints; j++) {
      current_iterate.original_evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, number_variables + j);
      });
      matrix.finalize(number_variables + j);
   }
   DEBUG << "KKT matrix for least-square multipliers:\n" << matrix << '\n';

   /********************************/
   /* generate the right-hand side */
   /********************************/
   initialize_vector(rhs, 0.);
   // objective gradient
   current_iterate.original_evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += model.objective_sign * derivative;
   });
   // variable bound constraints
   for (size_t i = 0; i < number_variables; i++) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i] + current_iterate.multipliers.upper_bounds[i];
   }
   DEBUG << "RHS for least-square multipliers: "; print_vector(DEBUG, rhs, 0, matrix.dimension);

   /********************/
   /* solve the system */
   /********************/
   std::vector<double> solution(matrix.dimension);
   solver.factorize(matrix);
   solver.solve(matrix, rhs, solution);
   DEBUG << "Solution: "; print_vector(DEBUG, solution, 0, matrix.dimension); DEBUG << '\n';

   // if least-square multipliers too big, discard them. Otherwise, store them
   if (norm_inf(solution, Range(number_variables, number_variables + model.number_constraints)) <= multipliers_max_norm) {
      for (size_t j = 0; j < model.number_constraints; j++) {
         multipliers[j] = solution[number_variables + j];
      }
   }
}