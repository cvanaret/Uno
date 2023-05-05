// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Preprocessing.hpp"
#include "solvers/QP/BQPDSolver.hpp"
#include "linear_algebra/CSCSymmetricMatrix.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "tools/Range.hpp"

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
void Preprocessing::compute_least_square_multipliers(const Model& model, SymmetricMatrix<double>& matrix, std::vector<double>& rhs,
      SymmetricIndefiniteLinearSolver<double>& linear_solver, Iterate& current_iterate, std::vector<double>& multipliers, double multiplier_max_norm) {
   current_iterate.evaluate_objective_gradient(model);
   current_iterate.evaluate_constraint_jacobian(model);

   /* build the symmetric matrix */
   matrix.reset();
   // identity block
   for (size_t i: Range(model.number_variables)) {
      matrix.insert(1., i, i);
      matrix.finalize_column(i);
   }
   // Jacobian of general constraints
   for (size_t j: Range(model.number_constraints)) {
      current_iterate.evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, model.number_variables + j);
      });
      matrix.finalize_column(model.number_variables + j);
   }
   DEBUG << "Matrix for least-square multipliers:\n" << matrix << '\n';

   /* generate the right-hand side */
   initialize_vector(rhs, 0.);
   // objective gradient
   current_iterate.evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += model.objective_sign * derivative;
   });
   // variable bound constraints
   for (size_t i: Range(model.number_variables)) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i] + current_iterate.multipliers.upper_bounds[i];
   }
   DEBUG << "RHS for least-square multipliers: "; print_vector(DEBUG, rhs, 0, matrix.dimension);
   
   /* solve the system */
   std::vector<double> solution(matrix.dimension);
   linear_solver.factorize(matrix);
   linear_solver.solve_indefinite_system(matrix, rhs, solution);
   DEBUG << "Solution: "; print_vector(DEBUG, solution, 0, matrix.dimension);

   // if least-square multipliers too big, discard them. Otherwise, keep them
   if (norm_inf(solution, Range(model.number_variables, model.number_variables + model.number_constraints)) <= multiplier_max_norm) {
      for (size_t j: Range(model.number_constraints)) {
         multipliers[j] = solution[model.number_variables + j];
      }
   }
   else {
      DEBUG << "Ignoring the least-square multipliers\n";
   }
   DEBUG << '\n';
}

size_t count_infeasible_linear_constraints(const Model& model, const std::vector<double>& constraint_values) {
   size_t infeasible_linear_constraints = 0;
   for (size_t j: model.get_linear_constraints()) {
      if (constraint_values[j] < model.get_constraint_lower_bound(j) || model.get_constraint_upper_bound(j) < constraint_values[j]) {
         infeasible_linear_constraints++;
      }
   }
   INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";
   return infeasible_linear_constraints;
}

void Preprocessing::enforce_linear_constraints(const Options& options, const Model& model, std::vector<double>& x, Multipliers& multipliers) {
   const auto& linear_constraints = model.get_linear_constraints();
   INFO << "Preprocessing phase: the problem has " << linear_constraints.size() << " linear constraints\n";
   if (not linear_constraints.empty()) {
      // evaluate the constraints
      std::vector<double> constraints(model.number_constraints);
      model.evaluate_constraints(x, constraints);
      const size_t infeasible_linear_constraints = count_infeasible_linear_constraints(model, constraints);
      if (0 < infeasible_linear_constraints) {
         // Hessian
         const CSCSymmetricMatrix<double> hessian = CSCSymmetricMatrix<double>::identity(model.number_variables);
         // constraint Jacobian
         RectangularMatrix<double> constraint_jacobian(linear_constraints.size());
         for (auto& constraint_gradient: constraint_jacobian) {
            constraint_gradient.reserve(model.number_variables);
         }
         for (size_t linear_constraint_index: Range(linear_constraints.size())) {
            const size_t j = linear_constraints[linear_constraint_index];
            model.evaluate_constraint_gradient(x, j, constraint_jacobian[linear_constraint_index]);
         }
         // variables bounds
         std::vector<Interval> variables_bounds(model.number_variables);
         for (size_t i: Range(model.number_variables)) {
            variables_bounds[i] = {model.get_variable_lower_bound(i) - x[i], model.get_variable_upper_bound(i) - x[i]};
         }
         // constraints bounds
         std::vector<Interval> constraints_bounds(linear_constraints.size());
         for (size_t linear_constraint_index: Range(linear_constraints.size())) {
            const size_t j = linear_constraints[linear_constraint_index];
            constraints_bounds[linear_constraint_index] =
                  {model.get_constraint_lower_bound(j) - constraints[j], model.get_constraint_upper_bound(j) - constraints[j]};
         }

         // solve the strictly convex QP
         BQPDSolver solver(model.number_variables, linear_constraints.size(), model.number_variables, true, options);
         std::vector<double> d0(model.number_variables); // = 0
         SparseVector<double> linear_objective; // empty
         WarmstartInformation warmstart_information{true, true, true, true};
         Direction direction = solver.solve_QP(model.number_variables, linear_constraints.size(), variables_bounds, constraints_bounds,
               linear_objective, constraint_jacobian, hessian, d0, warmstart_information);
         if (direction.status == SubproblemStatus::INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         // take the step
         add_vectors(x, direction.primals, 1., x);
         add_vectors(multipliers.lower_bounds, direction.multipliers.lower_bounds, 1., multipliers.lower_bounds);
         add_vectors(multipliers.upper_bounds, direction.multipliers.upper_bounds, 1., multipliers.upper_bounds);
         for (size_t linear_constraint_index: Range(linear_constraints.size())) {
            const size_t j = linear_constraints[linear_constraint_index];
            multipliers.constraints[j] += direction.multipliers.constraints[linear_constraint_index];
         }
         DEBUG2 << "Linear feasible initial point: "; print_vector(DEBUG2, x);
      }
   }
}