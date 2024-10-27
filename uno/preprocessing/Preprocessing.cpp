// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Preprocessing.hpp"
#include "ingredients/subproblems/Direction.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/CSCSparseStorage.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "solvers/QPSolver.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   // compute a least-square approximation of the multipliers by solving a linear system
   void Preprocessing::compute_least_square_multipliers(const Model& model, SymmetricMatrix<size_t, double>& matrix, Vector<double>& rhs,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver, Iterate& current_iterate, Vector<double>& multipliers,
         double multiplier_max_norm) {
      current_iterate.evaluate_objective_gradient(model);
      current_iterate.evaluate_constraint_jacobian(model);
      DEBUG << "Computing least-square multipliers\n";
      DEBUG2 << "Current primals: " << current_iterate.primals << '\n';

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
      DEBUG2 << "RHS for least-square multipliers: "; print_vector(DEBUG2, view(rhs, 0, model.number_variables + model.number_constraints));

      // if the residuals on the RHS are all 0, the least-square multipliers are all 0
      if (norm_inf(view(rhs, 0, model.number_variables + model.number_constraints)) == 0.) {
         multipliers.fill(0.);
         DEBUG << "Least-square multipliers are all 0.\n";
         return;
      }

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

      /* solve the system */
      Vector<double> solution(matrix.dimension());
      linear_solver.factorize(matrix);
      linear_solver.solve_indefinite_system(matrix, rhs, solution);

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

   size_t count_infeasible_linear_constraints(const Model& model, const std::vector<double>& constraint_values) {
      size_t infeasible_linear_constraints = 0;
      for (size_t constraint_index: model.get_linear_constraints()) {
         if (constraint_values[constraint_index] < model.constraint_lower_bound(constraint_index) ||
             model.constraint_upper_bound(constraint_index) < constraint_values[constraint_index]) {
            infeasible_linear_constraints++;
         }
      }
      return infeasible_linear_constraints;
   }

   void Preprocessing::enforce_linear_constraints(const Model& model, Vector<double>& x, Multipliers& multipliers, QPSolver& qp_solver) {
      const auto& linear_constraints = model.get_linear_constraints();
      INFO << "Preprocessing phase: the problem has " << linear_constraints.size() << " linear constraints\n";
      if (not linear_constraints.empty()) {
         // evaluate the constraints
         std::vector<double> constraints(model.number_constraints);
         model.evaluate_constraints(x, constraints);
         const size_t infeasible_linear_constraints = count_infeasible_linear_constraints(model, constraints);
         INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";
         if (0 < infeasible_linear_constraints) {
            // Hessian
            SymmetricMatrix<size_t, double> hessian(model.number_variables, model.number_variables, false, "CSC");
            for (size_t row_index: Range(model.number_variables)) {
               hessian.insert(1., row_index, row_index);
               hessian.finalize_column(row_index);
            }
            // constraint Jacobian
            RectangularMatrix<double> constraint_jacobian(linear_constraints.size(), model.number_variables);
            for (size_t linear_constraint_index: Range(linear_constraints.size())) {
               const size_t constraint_index = linear_constraints[linear_constraint_index];
               model.evaluate_constraint_gradient(x, constraint_index, constraint_jacobian[linear_constraint_index]);
            }
            // variable bounds
            std::vector<double> variables_lower_bounds(model.number_variables);
            std::vector<double> variables_upper_bounds(model.number_variables);
            for (size_t variable_index: Range(model.number_variables)) {
               variables_lower_bounds[variable_index] = model.variable_lower_bound(variable_index) - x[variable_index];
               variables_upper_bounds[variable_index] = model.variable_upper_bound(variable_index) - x[variable_index];
            }
            // constraint bounds
            std::vector<double> constraints_lower_bounds(linear_constraints.size());
            std::vector<double> constraints_upper_bounds(linear_constraints.size());
            for (size_t linear_constraint_index: Range(linear_constraints.size())) {
               const size_t constraint_index = linear_constraints[linear_constraint_index];
               constraints_lower_bounds[linear_constraint_index] = model.constraint_lower_bound(constraint_index) - constraints[constraint_index];
               constraints_upper_bounds[linear_constraint_index] = model.constraint_upper_bound(constraint_index) - constraints[constraint_index];
            }

            // solve the strictly convex QP
            Vector<double> d0(model.number_variables); // = 0
            SparseVector<double> linear_objective; // empty
            WarmstartInformation warmstart_information{true, true, true, true};
            Direction direction(model.number_variables, model.number_constraints);
            qp_solver.solve_QP(model.number_variables, linear_constraints.size(), variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
                  constraints_upper_bounds, linear_objective, constraint_jacobian, hessian, d0, direction, warmstart_information);
            if (direction.status == SubproblemStatus::INFEASIBLE) {
               throw std::runtime_error("Linear constraints cannot be satisfied at the initial point");
            }

            // take the step
            x += direction.primals;
            multipliers.lower_bounds += direction.multipliers.lower_bounds;
            multipliers.upper_bounds += direction.multipliers.upper_bounds;
            for (size_t linear_constraint_index: Range(linear_constraints.size())) {
               const size_t constraint_index = linear_constraints[linear_constraint_index];
               multipliers.constraints[constraint_index] += direction.multipliers.constraints[linear_constraint_index];
            }
            DEBUG3 << "Linear feasible initial point: " << x << '\n';
         }
      }
   }
} // namespace