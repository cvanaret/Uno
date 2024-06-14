// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <vector>
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "symbolic/VectorView.hpp"

// forward declarations
class Iterate;
class Model;
struct Multipliers;
class QPSolver;
template <typename IndexType, typename ElementType>
class SymmetricIndefiniteLinearSolver;
template <typename ElementType>
class Vector;

class Preprocessing {
public:
   template <typename SymmetricMatrix, typename LinearSolver>
   static void compute_least_square_multipliers(const Model& model, SymmetricMatrix& matrix, Vector<double>& rhs,
         LinearSolver& linear_solver, Iterate& current_iterate, Vector<double>& multipliers,
         double multiplier_max_norm);
   [[nodiscard]] static bool enforce_linear_constraints(const Model& model, Vector<double>& x, Multipliers& multipliers, QPSolver& qp_solver);
};

// compute a least-square approximation of the multipliers by solving a linear system
template <typename SymmetricMatrix, typename LinearSolver>
void Preprocessing::compute_least_square_multipliers(const Model& model, SymmetricMatrix& matrix, Vector<double>& rhs,
      LinearSolver& linear_solver, Iterate& current_iterate, Vector<double>& multipliers, double multiplier_max_norm) {
   current_iterate.evaluate_objective_gradient(model);
   current_iterate.evaluate_constraint_jacobian(model);

   /* build the symmetric matrix */
   matrix.reset();
   // identity block
   for (size_t variable_index: Range(model.number_variables)) {
      matrix.insert(1., variable_index, variable_index);
      matrix.finalize_column(variable_index);
   }
   // Jacobian of general constraints
   for (size_t constraint_index: Range(model.number_constraints)) {
      for (const auto [variable_index, derivative]: current_iterate.evaluations.constraint_jacobian[constraint_index]) {
         matrix.insert(derivative, variable_index, model.number_variables + constraint_index);
      }
      matrix.finalize_column(model.number_variables + constraint_index);
   }
   DEBUG2 << "Matrix for least-square multipliers:\n" << matrix << '\n';

   /* generate the right-hand side */
   rhs.fill(0.);
   // objective gradient
   for (const auto [variable_index, derivative]: current_iterate.evaluations.objective_gradient) {
      rhs[variable_index] += model.objective_sign * derivative;
   }
   // variable bound constraints
   for (size_t variable_index: Range(model.number_variables)) {
      rhs[variable_index] -= current_iterate.multipliers.lower_bounds[variable_index] + current_iterate.multipliers.upper_bounds[variable_index];
   }
   DEBUG2 << "RHS for least-square multipliers: "; print_vector(DEBUG2, view(rhs, 0, matrix.dimension));

   /* solve the system */
   Vector<double> solution(matrix.dimension);
   linear_solver.solve_indefinite_system(matrix, rhs, solution, true);

   // if least-square multipliers too big, discard them. Otherwise, keep them
   const auto trial_multipliers = view(solution, model.number_variables, model.number_variables + model.number_constraints);
   DEBUG2 << "Trial multipliers: "; print_vector(DEBUG2, trial_multipliers);
   if (norm_inf(trial_multipliers) <= multiplier_max_norm) {
      multipliers = trial_multipliers;
   }
   else {
      DEBUG << "Ignoring the least-square multipliers\n";
   }
   DEBUG << '\n';
}

#endif //UNO_PREPROCESSING_H
